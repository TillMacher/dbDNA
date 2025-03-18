import os
import sys
import re
import ast
import math
import glob
import gzip
import time
import shutil
import pickle
import threading
import random
import multiprocessing
import subprocess
from datetime import datetime
from pathlib import Path
import xmltodict
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from tqdm import tqdm
from joblib import Parallel, delayed
from requests_html import HTMLSession
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
import requests
import PyPDF2
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ete3 import NCBITaxa, Tree
from playwright.sync_api import sync_playwright

print_lock = threading.Lock()

####### STORE RUNTIMES #######

# Record the start time
t0 = datetime.now()

def time_diff(t0):
    t1 = datetime.now()
    # Calculate the elapsed time
    elapsed_time = t1 - t0
    # Format the elapsed time for readability
    hours, remainder = divmod(elapsed_time.seconds, 3600)
    minutes, seconds = divmod(remainder, 60)

    return f"{hours}h {minutes}min {seconds}s"

####### NEW PROJECT #######

def create_new_project(project_name, output_folder):

    # Main directory
    output_directory = Path('{}/{}_BarCodeBank'.format(output_folder, project_name))
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
        print('{} - Project "{}" was created.\n'.format(datetime.now().strftime("%H:%M:%S"), project_name))
    else:
        print('{} - Project "{}" already exists.\n'.format(datetime.now().strftime("%H:%M:%S"), project_name))

    # Sub directories for output and temporary files
    dir_1 = Path('{}/{}'.format(output_directory, '1_records_download'))
    if not os.path.exists(dir_1):
        os.makedirs(dir_1)

    dir_2 = Path('{}/{}'.format(output_directory, '2_phylogeny'))
    if not os.path.exists(dir_2):
        os.makedirs(dir_2)

    dir_3 = Path('{}/{}'.format(output_directory, '3_BarCodeBank'))
    if not os.path.exists(dir_3):
        os.makedirs(dir_3)

    return [output_directory, dir_1, dir_2, dir_3]

####### EXTRACT #######

## function to split the raw barcode table into families
def split_raw_barcode_table(output_directories, df_filtered):
    families = sorted(df_filtered['family'].drop_duplicates().values.tolist())

    for family in families:
        # Create subdirectories for each family
        family_dir = Path(f"{output_directories[2]}/{family}")
        family_dir.mkdir(parents=True, exist_ok=True)

        # Filter and ensure columns are string to avoid conversion issues
        sub_df = df_filtered[df_filtered['family'] == family].astype(str)

        # Save as parquet
        family_table = family_dir / f"{family}_raw_barcodes.parquet.snappy"
        sub_df.to_parquet(family_table, compression='snappy')

        print(f"{datetime.now():%H:%M:%S} - Created raw barcode table for {family}.")

## BOLD

def boldsystems_api(taxon, dir_out):
    taxon = ' '.join(taxon.split(' ')[:2])
    query_taxon = taxon.replace(" ", "%20")
    output = Path(f'{dir_out}/{taxon}.json')

    if Path(f'{output}.gz').is_file():
        tqdm.write(f'{datetime.now().strftime("%H:%M:%S")} - File already exists: {taxon}')
        return

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        context = browser.new_context(accept_downloads=True)
        page = context.new_page()

        url = f'https://portal.boldsystems.org/taxon/{query_taxon}'
        tqdm.write(f'{datetime.now().strftime("%H:%M:%S")} - Downloading from {url}')

        try:
            response = page.goto(url, wait_until="domcontentloaded")
            if not response or response.status != 200:
                tqdm.write(f'{datetime.now().strftime("%H:%M:%S")} - No BOLD records found for {taxon}')
                return

            download_button = page.query_selector('a[href*="json"], button:has-text("JSON")')
            if download_button:
                # Use expect_download to wait for download completion
                with page.expect_download() as download_info:
                    download_button.click()

                download = download_info.value  # Retrieve download object after completion
                download.save_as(output)
            else:
                tqdm.write(f'{datetime.now().strftime("%H:%M:%S")} - Download button not found.')

        except Exception as e:
            tqdm.write(f'{datetime.now().strftime("%H:%M:%S")} - An error occurred: {e}')

        finally:
            browser.close()
            if output.is_file():
                with open(output, 'rb') as f_in:
                    with gzip.open(f'{output}.gz', 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                os.remove(output)
                tqdm.write(f'{datetime.now().strftime("%H:%M:%S")} - Finished download for: {taxon}')

def download_data_from_bold(taxa_list, output_directories, marker):

    print('{} - Starting data download from bold.'.format(datetime.now().strftime("%H:%M:%S")))

    taxa_list_df = pd.read_excel(taxa_list).fillna('')
    dir_out = output_directories[1]
    taxa = sorted(set([i[0] for i in taxa_list_df.values.tolist() if i[0] != '']))
    max_retries = 3

    # Loop through taxa list and retry failed downloads if necessary
    for taxon in tqdm(taxa ,desc='Download', leave=True):
        retries = 0
        while retries < max_retries:
            try:
                boldsystems_api(taxon, dir_out)
                break  # Exit retry loop on success
            except Exception as e:
                retries += 1
                tqdm.write(f'{datetime.now().strftime("%H:%M:%S")} - Retry {retries}/{max_retries} for {taxon}: {e}')
                time.sleep(2)  # Short delay between retries

    print('{} - Finished data download from bold.\n'.format(datetime.now().strftime("%H:%M:%S")))

def extract_bold_json(output_directories, marker):

    print('{} - Starting data extraction from bold files.'.format(datetime.now().strftime("%H:%M:%S")))

    folder = output_directories[1]
    all_records = []

    for file in glob.glob(f'{folder}/*.json.gz'):
        df = pd.read_json(file, compression='gzip', lines=True).fillna('')

        # Ensure required columns are present, fill missing with empty strings if needed
        required_cols = {'species', 'family', 'nuc'}
        missing_cols = required_cols - set(df.columns)
        for col in missing_cols:
            df[col] = ''

        # Filter rows
        df = df[(df['species'] != '') & (df['family'] != '') & (df['nuc'] != '')]

        if not df.empty:
            all_records.append(df)
        else:
            print(f"{datetime.now():%H:%M:%S} - Warning: No species data found for {Path(file).stem} - species is removed.")

    # Concatenate all records into a single DataFrame, align columns as necessary
    df_filtered = pd.concat(all_records, ignore_index=True)
    df_filtered = df_filtered.loc[df_filtered['marker_code'] == marker]

    ## split and save a separate table for each family (will reduce runtimes significantly
    split_raw_barcode_table(output_directories, df_filtered)

    print('{} - Checkpoint: {}'.format(datetime.now().strftime("%H:%M:%S"), time_diff(t0)))
    print('{} - Finished data extraction from bold files.\n'.format(datetime.now().strftime("%H:%M:%S")))

## NCBI

def get_desired_ranks(taxid, desired_ranks):
    ncbi = NCBITaxa()
    lineage = ncbi.get_lineage(taxid)
    lineage2ranks = ncbi.get_rank(lineage)
    ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
    return {'{}_id'.format(rank): ranks2lineage.get(rank, '<not present>') for rank in desired_ranks}

def ncbi_taxid_request(taxid):

    desired_ranks = ['phylum', 'class', 'order', 'family', 'genus', 'species']
    taxonomy_list = []
    try:
        results = get_desired_ranks(taxid, desired_ranks)
        taxids = [str(taxid) for taxid in list(results.values())]

        # if the taxonomy is not present
        # DO THIS
        if '<not present>' in taxids:
            for taxid in taxids:
                if taxid != '<not present>':
                    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&id=' + str(taxid)
                    response = requests.get(url)
                    data = xmltodict.parse(response.content)
                    for entry in data['eSummaryResult']['DocSum']['Item']:
                        if entry['@Name'] == 'ScientificName':
                            name = entry['#text']
                            taxonomy_list.append(name)
                    time.sleep(0.2)
                else:
                    taxonomy_list.append('')

        # if all taxonomy information is present
        # DO THIS
        else:
            url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&id=' + ','.join(taxids)
            response = requests.get(url)
            data = xmltodict.parse(response.content)
            for entry in data['eSummaryResult']['DocSum']:
                for item in entry['Item']:
                    if item['@Name'] == 'ScientificName':
                        name = item['#text']
                        taxonomy_list.append(name)
        return taxonomy_list
    except ValueError:
        return ['No Match'] * 6

def accession2taxid(accession):
    url = 'https://www.ncbi.nlm.nih.gov/nuccore/{}'.format(accession)
    as_session = HTMLSession()
    as_session.headers.update({'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/98.9.4758.82 Safari/537.36'})
    retry_strategy = Retry(total = 10, status_forcelist = [400, 401, 403, 404, 429, 500, 502, 503, 504], backoff_factor = 1)
    adapter = HTTPAdapter(max_retries = retry_strategy)
    as_session.mount('https://', adapter)
    as_session.mount('http://', adapter)
    r = as_session.get(url, timeout = 300)
    data = r.text.split(';')
    taxid = [i for i in data if '?ORGANISM' in i][0].split('?')[-1].replace('&amp', '').replace('ORGANISM=', '')

    return taxid

def extract_MIDORI2_file(midori2_fasta, output_directories, taxa_list):

    print('{} - Starting data download from GenBank.\n'.format(datetime.now().strftime("%H:%M:%S")))

    download_folder = output_directories[1]
    taxa_list_all = pd.read_excel(taxa_list)['Taxon'].values.tolist()

    n_sequences = 0
    for record in SeqIO.parse(midori2_fasta, "fasta"):
        n_sequences += 1

    count = 0
    with open(midori2_fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            count += 1
            header = record.id
            accession = header.split(';')[0].split('.')[0]

            # Write the data to a .gb file
            gb_file = f"{download_folder}/{accession}.gb"

            if os.path.isfile(gb_file):
                print('{} - {} already exists ({}/{}).'.format(datetime.now().strftime("%H:%M:%S"), accession, count, n_sequences))

            elif any(species in header for species in taxa_list_all):

                # Always tell NCBI who you are
                Entrez.email = "till-hendrik.macher@uni-due.de"

                # Use Entrez.efetch to get the genbank record for the accession number
                handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")

                # Read the data returned
                data = handle.read()

                # write gb file
                with open(gb_file, "w") as f:
                    f.write(data)

                # Close the handle
                handle.close()

                # allow for maximum 3 requests per second
                time.sleep(1/3)

                print('{} - Finished {} ({}/{}).'.format(datetime.now().strftime("%H:%M:%S"), accession, count, n_sequences))

    print('{} - Finished data download from GenBank.\n'.format(datetime.now().strftime("%H:%M:%S")))

def extract_taxid(file):
    for record in SeqIO.parse(file, "genbank"):
        for feature in record.features:
            if feature.type == "source":
                taxid = feature.qualifiers.get("db_xref")
                if taxid:
                    for id in taxid:
                        if "taxon" in id:
                            return id.split(":")[1]

def extract_genbank_files(output_directories):

    print('{} - Starting to collect data from .gb files.\n'.format(datetime.now().strftime("%H:%M:%S")))

    taxids_xlsx = '/Volumes/Coruscant/dbDNA/taxids.xlsx'
    taxids_df = pd.read_excel(taxids_xlsx).fillna('No Match').drop_duplicates()
    taxids_dict = {i[0]:i[1:] for i in taxids_df.values.tolist()}

    files = glob.glob(f'{output_directories[1]}/*.gb')

    columns = ['processid', 'sampleid', 'recordID', 'catalognum', 'fieldnum', 'institution_storing', 'collection_code',
               'bin_uri', 'phylum_taxID', 'phylum_name', 'class_taxID', 'class_name', 'order_taxID', 'order_name', 'family_taxID',
               'family_name', 'subfamily_taxID', 'subfamily_name', 'genus_taxID', 'genus_name', 'species_taxID', 'species_name',
               'subspecies_taxID', 'subspecies_name', 'identification_provided_by', 'identification_method', 'identification_reference',
               'tax_note', 'voucher_status', 'tissue_type', 'collection_event_id', 'collectors', 'collectiondate_start',
               'collectiondate_end', 'collectiontime', 'collection_note', 'site_code', 'sampling_protocol', 'lifestage',
               'sex', 'reproduction', 'habitat', 'associated_specimens', 'associated_taxa', 'extrainfo', 'notes', 'lat',
               'lon', 'coord_source', 'coord_accuracy', 'elev', 'depth', 'elev_accuracy', 'depth_accuracy', 'country', 'province_state',
               'region', 'sector', 'exactsite', 'image_ids', 'image_urls', 'media_descriptors', 'captions', 'copyright_holders',
               'copyright_years', 'copyright_licenses', 'copyright_institutions', 'photographers', 'sequenceID', 'markercode',
               'genbank_accession', 'nucleotides', 'trace_ids', 'trace_names', 'trace_links', 'run_dates', 'sequencing_centers',
               'directions', 'seq_primers', 'marker_codes']

    all_references_list = []
    accession_taxid_taxonomy_list = []
    n_files = len(files)

    for i, file in enumerate(files):

        name = Path(file).stem

        reference_dict = {i: '' for i in columns}

        # open the gb file using SeqIO to scrape the basic information about each record
        for record in SeqIO.parse(file, "genbank"):
            reference_dict['nucleotides'] = str(record.seq)
            accession = record.id

            for feature in record.features:
                if feature.type == "source":
                    taxid = feature.qualifiers.get("db_xref")
                    if taxid:
                        for id in taxid:
                            if "taxon" in id:
                                taxid = int(id.split(":")[1])

            if taxid in taxids_dict.keys():
                taxonomy = taxids_dict[taxid]
                accession_taxid_taxonomy_list.append([accession] + taxonomy)
            else:
                # convert taxid to taxonomy
                taxonomy = ncbi_taxid_request(taxid)
                accession_taxid_taxonomy_list.append([accession] + taxonomy)
                taxids_dict[taxid] = taxonomy

        ## add data to dataframe
        reference_dict['phylum_name'] = taxonomy[0]
        reference_dict['class_name'] = taxonomy[1]
        reference_dict['order_name'] = taxonomy[2]
        reference_dict['family_name'] = taxonomy[3]
        reference_dict['genus_name'] = taxonomy[4]
        reference_dict['species_name'] = taxonomy[5]
        reference_dict['processid'] = f"Midori2-{i}"
        sampleid = Path(file).stem
        reference_dict['genbank_accession'] = sampleid
        reference_dict['sampleid'] = sampleid
        reference_dict['sequenceID'] = sampleid
        reference_dict['institution_storing'] = 'Mined from GenBank, NCBI'
        reference_dict['markercode'] = 'srRNA'

        # open the file again to scrape the remaining information
        # this is easier done line by line due to the missing information in many gb files
        with open(file, 'r') as f:
            for line in f:
                # country
                if 'country' in line:
                    try:
                        res = line.lstrip().rstrip('\n').replace('"', '').split('=')[1]
                        country = res.split(':')[0]
                        region = res.split(':')[1]
                        reference_dict['country'] = country
                        reference_dict['region'] = region
                    except:
                        reference_dict['country'] = line.lstrip().rstrip('\n').replace('"', '').replace("/country=", "")

                if 'AUTHORS' in line:
                    reference_dict['collectors'] = line.lstrip().rstrip('\n').replace('"', '').replace("AUTHORS", "").lstrip()

                # identifier
                if 'identified' in line:
                    reference_dict['identification_provided_by'] = line.lstrip().rstrip('\n').replace("/identified_by=", "").replace('\"', '')

                # lat lon
                if '/lat_lon' in line:
                    reference_dict['lat'] = line.lstrip().rstrip('\n').replace('/lat_lon=', '').replace('\"', '')

        print('{} - Finished {} ({}/{}).'.format(datetime.now().strftime("%H:%M:%S"), name, i+1, n_files))

        ## add reference sequences to list
        all_references_list.append(list(reference_dict.values()))

    ## create a dataframe
    df_filtered = pd.DataFrame(all_references_list, columns=columns)
    ## split and save a separate table for each family (will reduce runtimes significantly
    split_raw_barcode_table(output_directories, df_filtered)

    # update taxids file
    taxids_dict_2 = pd.DataFrame([[key] + values for key,values in taxids_dict.items()], columns=taxids_df.columns)
    taxids_dict_2.to_excel(taxids_xlsx, index=False)

    print('{} - Checkpoint: {}'.format(datetime.now().strftime("%H:%M:%S"), time_diff(t0)))
    print('{} - Finished to collect data from .gb files.\n'.format(datetime.now().strftime("%H:%M:%S")))

####### PHYLOGENY #######

# run alignment calculation
def run_mafft(family, mafft_executable, folder, outgroup_fasta, cpu_count, cpu_mode):

    # family = 'Philopotamidae'
    # outgroup_fasta = '/Volumes/Coruscant/dbDNA/outgroup.fasta'

    # Sub directories for output and temporary files
    family_dir = Path('{}/{}'.format(folder, family))
    if not os.path.exists(family_dir):
        os.makedirs(family_dir)

    ## required files to create
    fasta_file = Path('{}/{}.fasta'.format(family_dir, family))
    aln_file = Path('{}/{}.aln'.format(family_dir, family))
    species_file = Path('{}/{}.species'.format(family_dir, family))
    species_file_txt = Path(str(species_file) + ".txt")
    species_file_txt_snappy = Path(str(species_file_txt) + ".parquet.snappy")
    family_table = Path('{}/{}_raw_barcodes.parquet.snappy'.format(str(family_dir), family))
    ambiguous_barcodes_table = Path('{}/{}_ambiguous_taxa.parquet.snappy'.format(str(family_dir), family))
    usable_barcodes_table = Path('{}/{}_usable_taxa.parquet.snappy'.format(str(family_dir), family))
    df = pd.read_parquet(family_table)

    f_out = open(f'{family_dir}/0_mafft.stdout.txt', 'w')
    f_err = open(f'{family_dir}/0_mafft.stderr.txt', 'w')

    # Check if file already exists
    if os.path.isfile(species_file_txt_snappy):
        message = '{} - {} was already analysed!'.format(datetime.now().strftime("%H:%M:%S"),family)

    else:
        ## Create two files:
        ## ambiguous taxa (with ambiuous assignments)
        ## usable taxa (proper identification)
        raw_records = df.values.tolist()

        ## remove records with sp. and other special characters (i.e., ambiguous taxa)
        usable_records = []
        ambiguous_assignments = []
        for record in raw_records:
            # remove records with ambigous assignments
            species = record[21]
            special_characters = '!@#$%^&*()-+?_=,.<>/\'\"0123456789'
            # search for special characters in the species name
            if any(c in special_characters for c in species):
                ambiguous_assignments.append([record[4], record[21], '', '', 'Ambiguous assignment', -10, record[-19]])
            else:
                usable_records.append([record[4], record[21], '', '', '', 0, record[-19]])

        ## write ambiguous taxa to file
        ambiguous_table_df = pd.DataFrame(ambiguous_assignments, columns=['Sequence ID', 'Species Name', 'Cluster', 'Clade', 'State', 'Rating', 'Sequence'])
        ambiguous_table_df.to_parquet(ambiguous_barcodes_table)

        ## write ambiguous taxa to file
        usable_table_df = pd.DataFrame(usable_records, columns=['Sequence ID', 'Species Name', 'Cluster', 'Clade', 'State', 'Rating', 'Sequence'])
        usable_table_df.to_parquet(usable_barcodes_table)

        ## test if mafft must be run
        if os.path.isfile(aln_file) and os.path.getsize(aln_file) > 0:
            message = '{} - Alignment was already calculated for {}.'.format(datetime.now().strftime("%H:%M:%S"),family)

        # at least three sequences are required for the alignment!
        elif len(usable_records) < 4:
            message =  '{} - Skipped alignment for {} (less than 4 sequences).'.format(datetime.now().strftime("%H:%M:%S"),family)
            # create dummy file for all records
            f = open(species_file_txt, 'w')
            for record in raw_records:
                f.write('{}__{}\n'.format(record[68], record[21]))
            f.close()

        else:
            ############################################################################################################
            ## 1 - create fasta file

            # save outgroup names
            outgroup = []

            # create fasta file
            f = open(fasta_file, 'w')

            # enumerate to create unique identifier, which BOLD does not manage to do
            for i, record in enumerate(usable_records):
                ## format the fasta header
                header = '>{}__{}__{}\n'.format(record[0], record[1].replace(" ", '_'), i)
                ## write header
                f.write(header)
                ## write sequence
                sequence = '{}\n'.format(record[-1]).replace('-', 'N')
                f.write(sequence)

            ## also add the sequence from the outgroup file
            for i, record in enumerate(SeqIO.parse(outgroup_fasta, "fasta")):
                name = re.sub('\W+', '', record.id)
                header = f">{name}_outgroup_{i}\n"
                ## write header
                f.write(header)
                outgroup.append(f"{name}_outgroup_{i}")
                ## write sequence
                seq = f"{str(record.seq)}\n".replace('-', 'N')
                f.write(seq)

            f.close()

            ############################################################################################################
            ## 2 - create alignemnt using mafft

            # create alignment
            if cpu_mode == 'Parallel':
                command = "{} --auto --preservecase {} > {}".format(mafft_executable, fasta_file, aln_file)
            else:
                command = "{} --auto --thread {} --preservecase {} > {}".format(mafft_executable, cpu_count, fasta_file, aln_file)
            # Execute the command and capture stdout and stderr
            process = subprocess.Popen(command, shell=True, stdout=f_out, stderr=f_err)
            process.wait()

            message = '{} - Finished mafft alignment for {}.'.format(datetime.now().strftime("%H:%M:%S"), family)

    f_err.close()
    f_out.close()

    print(message)

## run phylogenetic tree calculation
def run_phylo(family, iqtree_executable, folder, cpu_count, cpu_mode):

    # family = 'Acrididae'
    # outgroup_fasta = '/Volumes/Coruscant/dbDNA/outgroup.fasta'

    # Sub directories for output and temporary files
    family_dir = Path('{}/{}'.format(folder, family))
    if not os.path.exists(family_dir):
        os.makedirs(family_dir)

    ## required files to create
    aln_file = Path('{}/{}.aln'.format(family_dir, family))
    tree_file = Path(str(aln_file) + '.treefile')
    species_file = Path('{}/{}.species'.format(family_dir, family))
    species_file_txt = Path(str(species_file) + ".txt")
    species_file_txt_snappy = Path(str(species_file_txt) + ".parquet.snappy")
    family_table = Path('{}/{}_raw_barcodes.parquet.snappy'.format(str(family_dir), family))
    df = pd.read_parquet(family_table)

    f_out = open(f'{family_dir}/1_phylo.stdout.txt', 'w')
    f_err = open(f'{family_dir}/1_phylo.stderr.txt', 'w')

    # Check if file already exists
    if os.path.isfile(species_file_txt_snappy):
        message = '{} - {} was already analysed!'.format(datetime.now().strftime("%H:%M:%S"), family)

    elif not os.path.isfile(aln_file):
        message = '{} - {} is missing an alignment file. Cannot calculate tree!'.format(datetime.now().strftime("%H:%M:%S"), family)
        # create dummy species txt file
        raw_records = df.values.tolist()
        f = open(species_file_txt, 'w')
        for record in raw_records:
            f.write('{}__{}\n'.format(record[68], record[21]))
        f.close()

    else:
        ## create tree using iqtree2 fast
        if os.path.isfile(tree_file):
            message = '{} - {} already has a tree file!'.format(datetime.now().strftime("%H:%M:%S"), family)
        else:
            ## calculate tree using iqtree2
            if cpu_mode == 'Parallel':
                command = [iqtree_executable, '-s', str(aln_file), '-m', 'K2P', '--ninit', '1', '--fast', '--redo', '-t', 'PARS']
            else:
                command = [iqtree_executable, '-s', str(aln_file), '-m', 'K2P', '-T', 'AUTO', '--threads-max', str(cpu_count), '--ninit', '1', '--fast', '--redo', '-t', 'PARS']
            subprocess.call(command, stdout=f_out, stderr=f_err)

            message = '{} - Finished phylogenetic tree for {}.'.format(datetime.now().strftime("%H:%M:%S"),family)

    f_err.close()
    f_out.close()

    print(message)

## run species delimitation
def run_mptp(family, mptp_executable, folder, outgroup_fasta):

    # family = 'Haliplidae'

    # Sub directories for output and temporary files
    family_dir = Path('{}/{}'.format(folder, family))
    if not os.path.exists(family_dir):
        os.makedirs(family_dir)

    ## required files to create
    aln_file = Path('{}/{}.aln'.format(family_dir, family))
    tree_file = Path(str(aln_file) + '.treefile')
    mldist_file = Path(str(aln_file) + '.mldist')
    minbr_file = Path(str(aln_file) + '.minbr')
    iqtree_log = Path(str(aln_file) + '.log')
    species_file = Path('{}/{}.species'.format(family_dir, family))
    species_file_txt = Path(str(species_file) + ".txt")
    species_file_txt_snappy = Path(str(species_file_txt) + ".parquet.snappy")
    family_table = Path('{}/{}_raw_barcodes.parquet.snappy'.format(str(family_dir), family))
    ambiguous_barcodes_table = Path('{}/{}_ambiguous_taxa.parquet.snappy'.format(str(family_dir), family))
    df = pd.read_parquet(family_table)

    ## collect bold records
    raw_records = df.values.tolist()

    ## also add the sequence from the outgroup file
    outgroup = []
    for i, record in enumerate(SeqIO.parse(outgroup_fasta, "fasta")):
        name = re.sub('\W+', '', record.id)
        outgroup.append(f"{name}_outgroup_{i}")

    ## open log files
    f_out = open(f'{family_dir}/2_mptp.stdout.txt', 'w')
    f_err = open(f'{family_dir}/2_mptp.stderr.txt', 'w')

    # Check if file already exists
    if os.path.isfile(species_file_txt_snappy):
        message = '{} - {} was already analysed!'.format(datetime.now().strftime("%H:%M:%S"),family)

    elif not os.path.isfile(tree_file):
        message = '{} - {} is missing an alignment file. Cannot calculate tree!'.format(datetime.now().strftime("%H:%M:%S"), family)
        # create dummy species txt file
        f = open(species_file_txt, 'w')
        for record in raw_records:
            f.write('{}__{}\n'.format(record[4], record[21]))
        f.close()

    else:
        # Estimation of minimum branch length
        command = f"{mptp_executable} --tree_file {str(tree_file)} --minbr_auto {str(aln_file)} --output_file {str(minbr_file)}"
        subprocess.run(command, shell=True, stdout=f_out, stderr=f_err)

        # Regular expression to find the minbr value
        minbr_pattern = r"--minbr\)\s+should be set to ([0-9]*\.?[0-9]+)"
        # Read the file and extract the minbr value
        with open(str(minbr_file) + '.txt', 'r') as file:
            content = file.read()
            match = re.search(minbr_pattern, content)
            if match:
                minbr_value = float(match.group(1))  # Convert the extracted value to float
            else:
                minbr_value = 0

        # identify species using mptp
        command = f"{mptp_executable} --ml --multi --minbr {minbr_value} --tree_file {str(tree_file)} --outgroup {outgroup[0]} --outgroup_crop --output {str(species_file)}"
        subprocess.run(command, shell=True, stdout=f_out, stderr=f_err)

        ## sleep to give time to write the file (happened several times that the file was not yet properly written)
        time.sleep(10)

        ## check if a species delimitation file was created
        if not os.path.isfile(species_file_txt):
            # create dummy species txt file
            f = open(species_file_txt, 'w')
            for record in raw_records:
                f.write('{}__{}\n'.format(record[0], record[1]))
            f.close()
            message = '{} - mptp encountered an error for {}.'.format(datetime.now().strftime("%H:%M:%S"), family)

        else:
            ## MPTP was run successfully and the file can be converted!

            # collect data
            current_species = None
            species_data = []
            sp = 0

            # Open and read the file
            with open(species_file_txt, 'r') as file:
                for line in file:
                    if 'LRT: passed' in line:
                        LRT = 'Passed'
                    if 'LRT: failed' in line:
                        LRT = 'Failed'

                    # Skip the initial information
                    if line.startswith('Species') and line.endswith(':\n'):
                        sp += 1
                        current_species = '{}_{}'.format(family, sp)
                    elif current_species != None and line != '\n':
                        # Parse the species data
                        ## collect the temporary id
                        ids = line.strip().split(';;')

                        ## loop through all ids (can be more than one records, depending on the duplicates)
                        for record in ids:
                            ## collect original sequence information
                            species_id = record.split('__')[0]
                            species_name = record.split('__')[1]
                            species_data.append([species_id, species_name, current_species, LRT])

                ## also append the ambiguous taxa here
                ambiguous_taxa = pd.read_parquet(ambiguous_barcodes_table).values.tolist()
                for record in ambiguous_taxa:
                    species_id = record[0]
                    species_name = record[1]
                    current_species = ''
                    LRT = 'Ambiguous identification'
                    species_data.append([species_id, species_name, current_species, LRT])

            ## check if IQ-TREE worked properly... if not fall back to initial tree
            if not os.path.isfile(mldist_file):
                print('{} - WARNING! IQ-Tree failed for {} - species delimitation might be biased.'.format(datetime.now().strftime("%H:%M:%S"), family))
                ## collect removed records
                identical_storage = []
                f = open(iqtree_log, 'r')
                for line in f:
                    line = line.strip()
                    if line.startswith('NOTE:') and 'is ignored but added at the end' in line:
                        species_a = line.split()[1]
                        species_b = line.split()[4].rstrip(')')
                        # a is ignored (because it is identical to b) and must be re-added
                        identical_storage.append([species_a, species_b])
                f.close()
                identical_storage_df = pd.DataFrame(identical_storage, columns=['Identical', 'Kept'])

                ## collect information of identical records
                for removed, kept in identical_storage_df.values.tolist():
                    ## collect information about REMOVED identical record
                    species_id = removed.split('__')[0]
                    species_name = removed.split('__')[1]
                    ## collect information about KEPT identical record
                    ref_id = kept.split('__')[0]
                    ref_name = kept.split('__')[1]
                    res = [i for i in species_data if ref_id in i and ref_name in i][0]
                    if res != []:
                        current_species = res[-2]
                        LRT = res[-1]
                    else:
                        current_species = 'NA'
                        LRT = 'ERROR'
                    species_data.append([species_id, species_name, current_species, LRT])

            ## then create a species delimitation table
            species_df = pd.DataFrame(species_data, columns=['Sequence ID', 'Species Name', 'Cluster', 'LRT'])
            res = []
            for species_group in set(species_df['Cluster'].values.tolist()):
                sub_df = species_df.loc[species_df['Cluster'] == species_group]
                LRT = list(set(sub_df['LRT'].values.tolist()))[0]
                sub_df = sub_df.drop('LRT', axis=1)
                n_records = len(sub_df)
                n_species = len(set(sub_df['Species Name'].values.tolist()))

                ## monophyletic and more than 2 sequences == monophyletic (1)
                if LRT == 'Passed' and n_records >= 2 and n_species == 1:
                    phylogeny = 'monophyletic'

                ## monophyletic but single sequence == monophyletic (2)
                elif LRT == 'Passed' and n_records < 2 and n_species == 1:
                    phylogeny = 'monophyletic (singleton)'

                elif LRT == 'Passed' and n_species != 1:
                    phylogeny = 'paraphyletic'

                ## LRT failed == paraphyletic (1)
                elif LRT == 'Failed':
                    phylogeny = 'paraphyletic (LRT failed)'

                ## Mark ambiguous assignments
                elif LRT == 'Ambiguous identification':
                    phylogeny = 'Ambiguous identification'

                ## LRT insufficient data == paraphyletic (2)
                else:
                    phylogeny = 'paraphyletic (insufficient data)'

                ## add to results
                for record in sub_df.values.tolist():
                    res.append(record + [0, phylogeny, 0])

            ## write dataframe
            df_out = pd.DataFrame(res, columns=['Sequence ID', 'Species Name', 'Cluster', 'Clade', 'State', 'Rating'])
            df_out = df_out.astype('string')  # Convert to string, otherwise parquet crashes
            df_out.to_parquet(species_file_txt_snappy, index=False)
            message = '{} - Finished species delimitation for {}'.format(datetime.now().strftime("%H:%M:%S"), family)

    ####################################################################################################################
    ## create the species delimitation file

    ## check if a species delimitation file was created
    if not os.path.isfile(species_file_txt_snappy):
        ## write a dummy file if mptp failed
        species_data = [[record[68], record[21], '', '', 'Error', 0] for record in raw_records]
        species_df = pd.DataFrame(species_data, columns=['Sequence ID', 'Species Name', 'Cluster', 'Clade', 'State', 'Rating'])
        ## write dataframe
        species_df = species_df.astype('string')  # Convert to string, otherwise parquet crashes
        species_df.to_parquet(species_file_txt_snappy, index=False)
        message = '{} - Finished {}, but failed species delimitation'.format(datetime.now().strftime("%H:%M:%S"), family)

    f_out.close()
    f_err.close()

    print(message)

# main script for phylogenetic approach
def phylogenetic_approach(output_directories, mafft_executable, iqtree_executable, mptp_executable, cpu_count, cpu_mode):

    print('{} - Starting phylogenetic approach.'.format(datetime.now().strftime("%H:%M:%S")))

    # Store fasta, alignments and tree
    folder = output_directories[2]

    # Extract families
    files = glob.glob('{}/*/*_raw_barcodes.parquet.snappy'.format(output_directories[2]))
    families = sorted([Path(i).name.replace('_raw_barcodes.parquet.snappy', '') for i in files])

    print('{} - Using {}/{} CPUs.\n'.format(datetime.now().strftime("%H:%M:%S"), cpu_count, cpu_count+1))

    ####################################################################################################################
    # 1) create alignments using parallel and a single core per family
    if cpu_mode == 'Parallel':
        Parallel(n_jobs=cpu_count, backend='threading')(delayed(run_mafft)(family, mafft_executable, folder, outgroup_fasta, cpu_count, cpu_mode) for family in families)
    else:
        [run_mafft(family, mafft_executable, folder, outgroup_fasta, cpu_count, cpu_mode) for family in families]

    print('{} - Checkpoint: {}'.format(datetime.now().strftime("%H:%M:%S"), time_diff(t0)))
    print('{} - Finished mafft alignments.\n'.format(datetime.now().strftime("%H:%M:%S")))

    ####################################################################################################################
    # 2) calculate trees using all available cores per family - ALSO IN PARALLEL SINCE IQTREE CRASHES IN MUTLICORE MODE
    if cpu_mode == 'Parallel':
        Parallel(n_jobs = cpu_count, backend='threading')(delayed(run_phylo)(family, iqtree_executable, folder, cpu_count, cpu_mode) for family in families)
    else:
        [run_phylo(family, iqtree_executable, folder, cpu_count, cpu_mode) for family in families]

    print('{} - Checkpoint: {}'.format(datetime.now().strftime("%H:%M:%S"), time_diff(t0)))
    print('{} - Finished calculation of phylogenetic trees.\n'.format(datetime.now().strftime("%H:%M:%S")))

    ####################################################################################################################
    # 3) perform species delimitation
    Parallel(n_jobs = cpu_count, backend='threading')(delayed(run_mptp)(family, mptp_executable, folder, outgroup_fasta) for family in families)
    print('{} - Checkpoint: {}'.format(datetime.now().strftime("%H:%M:%S"), time_diff(t0)))
    print('{} - Finished species delimitation.\n'.format(datetime.now().strftime("%H:%M:%S")))

def haversine(lat1, lon1, lat2, lon2):
    # Convert latitude and longitude from degrees to radians
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])

    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
    c = 2 * math.asin(math.sqrt(a))

    # Radius of earth in kilometers is 6371
    km = 6371 * c
    return km

## collect all possiblitites for several metadata
def collect_metadata():
    files = glob.glob('/Volumes/Coruscant/dbDNA/FEI_genera_v2_BarCodeBank/1_records_download/*.tsv')
    test = ['identification_provided_by', 'identification_method', 'voucher_status']
    res = {}
    for file in files:
        df = pd.read_table(file, low_memory=False).fillna('')
        for t in test:
            if t not in res.keys():
                res[t] = [i for i in df[t].drop_duplicates().values.tolist() if i != '']
            else:
                res[t] = sorted(set(res[t] + [i for i in df[t].drop_duplicates().values.tolist() if i != '']))
    df1 = pd.DataFrame(res['identification_provided_by'], columns=['identification_provided_by'])
    df2 = pd.DataFrame(res['identification_method'], columns=['identification_method'])
    df3 = pd.DataFrame(res['voucher_status'], columns=['voucher_status'])

    ## filter method
    for value in df2['identification_method'].values.tolist():
        value = value.lower()
        keywords = ["bin", "bold", "tree"]
        exclude = ["morphology", "bold:"]
        if any(word in value for word in keywords) and not any(word in value for word in exclude):
            print(value)


    df1.to_excel('/Volumes/Coruscant/dbDNA/bold_info/identification_provided_by.xlsx', index=False)
    df2.to_excel('/Volumes/Coruscant/dbDNA/bold_info/identification_method.xlsx', index=False)
    df3.to_excel('/Volumes/Coruscant/dbDNA/bold_info/voucher_status.xlsx', index=False)

# function to rate each family
def rate_family(output_directories, identifier_whitelist_lst, location_whitelist_lst, family, use_coordinates, use_country, d1, d2, d3):

    # family = 'Acanthobdellidae'
    # check_coordinates = True
    # df  = pd.read_parquet('/Volumes/Coruscant/dbDNA/Small_db_BarCodeBank/2_phylogeny/Viviparidae/Viviparidae_raw_barcodes.parquet.snappy')
    # voucher_status
    # identification_method

    # Sub directories for output and temporary files
    family_dir = Path('{}/{}'.format(output_directories[2], family))
    if not os.path.exists(family_dir):
        os.makedirs(family_dir)

    ## relevant files
    species_file_txt_snappy = Path('{}/{}.species.txt.parquet.snappy'.format(str(family_dir), family))
    ratings_file = Path('{}/{}.ratings.parquet.snappy'.format(str(family_dir), family))
    family_table = Path('{}/{}_raw_barcodes.parquet.snappy'.format(str(family_dir), family))

    # Collect information about records
    df1 = pd.read_parquet(family_table).fillna('')

    ## store ratings in list
    all_ratings_list = []

    # Location white list
    main_country = [i for i in location_whitelist_lst['main country'].values.tolist() if i != '']
    neighbour_countries = [i for i in location_whitelist_lst['neighbour_countries'].values.tolist() if i != '']
    continent = [i for i in location_whitelist_lst['continent'].values.tolist() if i != '']

    ## read mPTP file
    df2 = pd.read_parquet(species_file_txt_snappy)

    ## rate each individual barcode of the family
    for record in df2.values.tolist():
        ## collect information about the record
        # family file
        record_id = record[0]
        species_group = record[2]
        phylogeny = record[4]
        rating = int(record[5])

        # raw table file
        raw_record = df1.loc[df1['record_id'] == record_id]
        record_id = raw_record['record_id'].values.tolist()[0]
        institution_storing = raw_record['inst'].values.tolist()[0]
        bin_uri = raw_record['bin_uri'].values.tolist()[0]
        phylum_name = raw_record['phylum'].values.tolist()[0]
        class_name = raw_record['class'].values.tolist()[0]
        order_name = raw_record['order'].values.tolist()[0]
        family_name = raw_record['family'].values.tolist()[0]
        genus_name = raw_record['genus'].values.tolist()[0]
        species_name = raw_record['species'].values.tolist()[0]
        identification_by = raw_record['identified_by'].values.tolist()[0]
        identification_method = raw_record['identification_method'].values.tolist()[0]
        reverse_bin = ''
        voucher_status = raw_record['voucher_type'].values.tolist()[0]
        country = raw_record['country/ocean'].values.tolist()[0]
        province = raw_record['province/state'].values.tolist()[0]
        region = raw_record['region'].values.tolist()[0]
        exactsite = raw_record['site'].values.tolist()[0]
        lifestage = raw_record['life_stage'].values.tolist()[0]
        sex = raw_record['sex'].values.tolist()[0]
        image_urls = ''
        markercode = raw_record['marker_code'].values.tolist()[0]
        nucleotides = raw_record['nuc'].values.tolist()[0]
        sequence_quality = 'Sufficient' ## default to "sufficient and only change when good or bad
        distance = ''

        ## rating_calc
        rating_calc = []

        ## phylogeny
        if phylogeny == 'monophyletic':
            rating += 15
            rating_calc.append('+15 (monophyletic)')
        elif phylogeny == 'monophyletic (singleton)':
            rating += 5
            rating_calc.append('+5 (monophyletic-singleton)')

        ## trimm leading and trailing gaps from barcode
        nucleotides_trimmed = nucleotides.strip('-').strip('N')
        barcode_length = len(nucleotides_trimmed)

        ## good sequence quality
        allowed_chars = set('ACGT')
        if set(nucleotides_trimmed).issubset(allowed_chars):
            rating += 5
            sequence_quality = 'Good'
            rating_calc.append('+5 (sequence quality)')

        ## bad sequence quality
        not_allowed_chars = [i for i in set(nucleotides_trimmed) if i not in allowed_chars]
        if len(not_allowed_chars) != 0:
            n_not_allowed_chars = sum([nucleotides_trimmed.count(i) for i in not_allowed_chars])
            rel = n_not_allowed_chars / len(nucleotides_trimmed) * 100
            if rel >= 2:
                rating -= 10
                sequence_quality = 'Bad'
                rating_calc.append('-10 (sequence quality)')

        ## reverse BIN taxonomy
        keywords = ["bin", "bold", "tree"]
        exclude = ["morphology", "bold:"]
        value = identification_method.lower()
        if any(word in value for word in keywords) and not any(word in value for word in exclude):
            rating -= 10
            reverse_bin = 'Investigate'
            rating_calc.append('-10 (rev. BIN identification)')

        ## sequence length (>= 500 bp are accepted as barcode)
        if len(nucleotides_trimmed) >= 500:
            rating += 5
            rating_calc.append('+5 (barcode length)')

        ## Identification white list
        if identification_by in identifier_whitelist_lst:
            rating += 10
            rating_calc.append('+10 (identifier on whiteilist)')

        if use_coordinates == 'yes':
            use_country = 'no'
            try:
                c1 = ast.literal_eval(raw_record['coord'].values.tolist()[0])[0]
                c2 = ast.literal_eval(raw_record['coord'].values.tolist()[0])[1]

                if c1 != '' and c2 == '':
                    lat = float(c1.split(' ')[0])
                    lon = float(c1.split(' ')[2])
                    if c1.split(' ')[-1] == 'W':
                        lon = -abs(lon)
                else:
                    lat = float(c1)
                    lon = float(c2)

                ## calculate distance to point of interest
                distance = haversine(lat, lon, lat_db, lon_db)
                if distance <= d1:
                    rating += 9
                    rating_calc.append('+9 (coordinates)')
                elif distance <= d2:
                    rating += 6
                    rating_calc.append('+6 (coordinates)')
                elif distance <= d3:
                    rating += 3
                    rating_calc.append('+3 (coordinates)')

            except ValueError:
                lat = raw_record['coord'].values.tolist()[0]
                lon = raw_record['coord'].values.tolist()[0]

        else:
            ## Sampling location with white list
            if country in main_country:
                rating += 9
                rating_calc.append('+9 (country)')
            elif country in neighbour_countries:
                rating += 6
                rating_calc.append('+6 (country)')
            elif country in continent:
                rating += 3
                rating_calc.append('+3 (country)')
            lat = raw_record['coord'].values.tolist()[0]
            lon = raw_record['coord'].values.tolist()[0]

        ## Check if photo is available
        if image_urls != '':
            rating += 1
            rating_calc.append('+1 (image)')

        ## Available metadata
        if province != '':
            rating += 1
            rating_calc.append('+1 (province)')
        if region != '':
            rating += 1
            rating_calc.append('+1 (region)')
        if exactsite != '':
            rating += 1
            rating_calc.append('+1 (exactsite)')
        if lifestage != '':
            rating += 1
            rating_calc.append('+1 (lifestage)')
        if sex != '':
            rating += 1
            rating_calc.append('+1 (sex)')

        ## always set to -20 if there is no proper identification
        if phylogeny == 'Ambiguous identification':
            rating = -20
            rating_calc.append('set -20 (Ambiguous identification)')

        rating_calc_str = '; '.join(rating_calc)

        all_ratings_list.append([rating, record_id, bin_uri, phylum_name, class_name, order_name,
                                 family_name, genus_name, species_name, phylogeny, species_group, identification_by, identification_method, reverse_bin, voucher_status,
                                 institution_storing, lat, lon, distance, country, province, region, exactsite, lifestage, sex, image_urls, markercode,
                                 sequence_quality, barcode_length, nucleotides, rating_calc_str])

    # create unfiltered dataframe
    ratings_df = pd.DataFrame(all_ratings_list, columns=["rating", "record_id", "bin_uri", "phylum_name", "class_name", "order_name",
                                 "family_name", "genus_name", "species_name", "phylogeny", "species_group", "identification_by", "identification_method", "reverse_bin", "voucher_status",
                                 "institution_storing", "lat", "lon", "distance", "country", "province", "region", "exactsite", "lifestage", "sex", "image_urls", "markercode",
                                 "sequence_quality", "barcode_length", "nucleotides", "rating calculation"])

    ratings_df = ratings_df.sort_values('rating', ascending=False)
    ratings_df = ratings_df.astype('string')  # Convert to string, otherwise parquet crashes
    ratings_df['rating'] = ratings_df['rating'].astype('int')
    #ratings_df['clade'] = ratings_df['clade'].astype('int')

    # write to parquet
    ratings_df.to_parquet(ratings_file, index=False)

    print('{} - Finished rating for {}.'.format(datetime.now().strftime("%H:%M:%S"), family))

# main script for rating algorithm
def rating_system(output_directories, identifier_whitelist, location_whitelist, project_name, cpu_count):

    print('{} - Collecting raw records.'.format(datetime.now().strftime("%H:%M:%S")))

    # Extract families
    files = glob.glob('{}/*/*_raw_barcodes.parquet.snappy'.format(output_directories[2]))
    families = sorted([Path(i).name.replace('_raw_barcodes.parquet.snappy', '') for i in files])

    identifier_whitelist_lst = sorted(set(pd.read_excel(identifier_whitelist, sheet_name='Identifier_Whitelist').fillna('')['Name (as used on BOLD)'].values.tolist()))
    location_whitelist_lst = pd.read_excel(location_whitelist).fillna('')

    ## for testing
    # files_to_delete = glob.glob('/Volumes/Coruscant/dbDNA/FEI_genera_BarCodeBank/2_phylogeny/*/*.ratings.parquet.snappy')
    # [os.remove(file) for file in files_to_delete if os.path.isfile(file)]

    print('{} - Starting to rate sequences.'.format(datetime.now().strftime("%H:%M:%S")))

    ## rate all families in parallel
    Parallel(n_jobs=cpu_count, backend='threading')(delayed(rate_family)(output_directories, identifier_whitelist_lst, location_whitelist_lst, family, use_coordinates, use_country, d1, d2, d3) for family in families)

    print('{} - Finished rating for all families.\n'.format(datetime.now().strftime("%H:%M:%S")))

    # Initialize an empty list to hold all dataframes
    dfs = []

    for family in families:
        family_dir = Path('{}/{}'.format(output_directories[2], family))
        ratings_snappy = family_dir / '{}.ratings.parquet.snappy'.format(family)

        if not ratings_snappy.is_file():
            print('{} - WARNING: No ratings file found for {}!!'.format(datetime.now().strftime("%H:%M:%S"), family))
        else:
            df = pd.read_parquet(ratings_snappy)
            dfs.append(df)  # Append the dataframe to the list

    # Concatenate all dataframes in the list
    ratings_df = pd.concat(dfs, ignore_index=True)

    ## sort by rating
    ratings_df = ratings_df.sort_values('rating', ascending=False)

    ## write all tables to a single file
    output_file_1 = Path('{}/{}.BarCodeBank.parquet.snappy'.format(output_directories[3], project_name))
    ratings_df.to_parquet(output_file_1, index=False)
    print('{} - Saved database to .parquet.snappy.'.format(datetime.now().strftime("%H:%M:%S")))

    ## test if the dataframe can be written to excel
    if ratings_df.shape[0] > 65000:
        output_file_1 = Path('{}/{}.BarCodeBank.csv'.format(output_directories[3], project_name))
        ratings_df.to_csv(output_file_1, index=False)
        print('{} - Unable to write to .xlsx. Saved database to .csv instead!'.format(datetime.now().strftime("%H:%M:%S")))
    else:
        output_file_1 = Path('{}/{}.BarCodeBank.xlsx'.format(output_directories[3], project_name))
        ratings_df.to_excel(output_file_1, index=False)
        print('{} - Saved database to .xlsx.'.format(datetime.now().strftime("%H:%M:%S")))

    ## write to .fasta file

    print('{} - Checkpoint: {}'.format(datetime.now().strftime("%H:%M:%S"), time_diff(t0)))
    print('{} - Finished to rate sequences.\n'.format(datetime.now().strftime("%H:%M:%S")))

# filter for specific target taxa (e.g., operational taxa list)
def filter_database(output_directories, filter_list):
    db_file = Path('{}/{}.BarCodeBank.parquet.snappy'.format(output_directories[3], project_name))
    ratings_df = pd.read_parquet(db_file)
    filter_list = '/Volumes/Coruscant/dbDNA/perlodes_TTT_conversion.xlsx'
    filter_df = pd.read_excel(filter_list).fillna('')

    OTL_df = pd.DataFrame()

    for taxon in filter_df.values.tolist():
        ## collect taxon
        query = [i for i in taxon[4:] if i != ''][0]
        for level in ['species_name', 'genus_name', 'family_name', 'order_name', 'class_name', 'phylum_name']:
            sub_df = ratings_df.loc[ratings_df[level] == query]
            if len(sub_df) != 0:
                break
        if len(sub_df) != 0:
            sub_df = sub_df.copy()
            sub_df.loc[:, 'OTL_taxon'] = query
            OTL_df = pd.concat([OTL_df, sub_df], ignore_index=True)


    ## removal of barcodes without identifier
    OTL_df_filtered_1 = OTL_df.loc[OTL_df['identification_by'] != '']

    ## removal of barcodes with potential reverse BIN taxonomy
    OTL_df_filtered_2 = OTL_df_filtered_1.loc[OTL_df_filtered_1['reverse_bin'] == '']

    ## drop OTL taxon
    OTL_df_filtered_3 = OTL_df_filtered_2.drop(['OTL_taxon'], axis=1)

    ## remove duplicates
    OTL_df_filtered_4 = OTL_df_filtered_3.drop_duplicates()

    len(OTL_df_filtered_4)

####### BLAST DATABASE #######

def create_database(output_directories, project_name, makeblastdb_exe):

    print('{}: Starting to build a new database.'.format(datetime.now().strftime('%H:%M:%S')))

    database_snappy = Path(f'{output_directories[3]}/{project_name}.BarCodeBank.parquet.snappy')
    df = pd.read_parquet(database_snappy)

    # Create fasta file
    fasta_file = Path(f'{output_directories[3]}/{project_name}.BarCodeBank.fasta')

    records = []
    for sequence in df.values.tolist():
        nucleotides = sequence[-1].replace('-', '')
        seq_id = sequence[1]
        process_id = sequence[2]
        species = sequence[9].replace(" ", "_")
        header = f'>{seq_id}__{process_id}__{species}'
        record = SeqRecord(Seq(nucleotides), id=header, description="")
        records.append(record)

    SeqIO.write(records, fasta_file, "fasta")

    # Create database
    ## collect files

    ## create a new folder for the database
    db_folder = Path('{}/4_{}_database'.format(output_directories[0], project_name))
    try:
        os.mkdir(db_folder)
    except FileExistsError:
        pass

    ## build a new database
    db_name = Path(db_folder).joinpath('db')
    subprocess.call([makeblastdb_exe, '-in', str(fasta_file), '-dbtype', 'nucl', '-out', str(db_name)])

    # Check if database was created
    print('{} - Checkpoint: {}'.format(datetime.now().strftime("%H:%M:%S"), time_diff(t0)))
    if os.path.isfile(f'{db_name}.ndb'):
        print('{}: Finished building database.'.format(datetime.now().strftime('%H:%M:%S')))
    else:
        print('{}: An error occurred when building the database.'.format(datetime.now().strftime('%H:%M:%S')))

####### ONLY FOR TESTING #######

def validate_database(output_directories):

    families = [Path(i).name for i in glob.glob(str(output_directories[2]) + "/*")]
    database_file = glob.glob(str(output_directories[3]) + "/*.BarCodeBank.parquet.snappy")[0]
    database_df = pd.read_parquet(database_file)
    final_processIDs = set(database_df['processid'].values.tolist())

    for family in tqdm(families):
        file = output_directories[2].joinpath(family, f'{family}_raw_barcodes.parquet.snappy')
        raw_table = pd.read_parquet(file)
        processIDs = set(raw_table['processid'].values.tolist())
        missing = len(processIDs - final_processIDs)
        if missing != 0:
            print(f'Missing records: {family}')

## delete all files in the phylogeny folder
def clear_folders():
    ## for testing
    files_to_delete = glob.glob('/Volumes/Coruscant/dbDNA/European_fish_BarCodeBank/2_phylogeny/*/*.ratings.parquet.snappy')
    [os.remove(file) for file in files_to_delete if os.path.isfile(file)]

    files_to_delete = glob.glob('/Volumes/Coruscant/dbDNA/FEI_genera_BarCodeBank/2_phylogeny/*/*.species.txt.parquet.snappy')
    [os.remove(file) for file in files_to_delete if os.path.isfile(file)]

########################################################################################################################

def optimal_height(y_values):
    n = len(y_values)
    min_n = 50
    min_height = 500
    if n > min_n:
        delta = n - min_n
        height = (delta * 20) + min_height
    else:
        height = min_height
    return height

def report(output_directories, project_name, taxa_list):
    ## import files
    barcode_bank_file = Path('{}/{}.BarCodeBank.parquet.snappy'.format(output_directories[3], project_name))
    reference_db_df = pd.read_parquet(barcode_bank_file).fillna('')

    # Extract families
    dfs = []
    files = glob.glob('{}/*/*_raw_barcodes.parquet.snappy'.format(output_directories[2]))
    families = sorted([Path(i).name.replace('_raw_barcodes.parquet.snappy', '') for i in files])
    for family in families:
        family_dir = Path('{}/{}'.format(output_directories[2], family))
        ratings_snappy = family_dir / '{}_raw_barcodes.parquet.snappy'.format(family)
        if not ratings_snappy.is_file():
            print('{} - WARNING: No ratings file found for {}!!'.format(datetime.now().strftime("%H:%M:%S"), family))
        else:
            df = pd.read_parquet(ratings_snappy)
            dfs.append(df)  # Append the dataframe to the list
    # Concatenate all dataframes in the list
    raw_data_df = pd.concat(dfs, ignore_index=True)

    ## check for missing reference barcodes
    ids_0 = raw_data_df['sequenceID'].drop_duplicates()
    ids_1 = reference_db_df['sequenceID'].drop_duplicates()
    len_0 = len(ids_0)
    len_1 = len(ids_1)
    n_shared = len(set(ids_0) & set(ids_1))
    shared_perc = n_shared / len_0 * 100

    ## 1) taxa without barcodes
    taxa_list_df = pd.read_excel(taxa_list)
    res = {}
    for taxon in taxa_list_df.values.tolist():
        taxon = taxon[0]
        res[taxon] = reference_db_df['genus_name'].values.tolist().count(taxon)
    barcodes_df = pd.DataFrame([[i,j] for i,j in res.items()], columns=['Taxon', 'Barcodes'])
    barcodes_df.to_excel(Path('{}/{}.report_barcodes.xlsx'.format(output_directories[3], project_name)), index=False)

    check_df = barcodes_df.copy()
    max_barcodes = max(barcodes_df['Barcodes'])
    max_rounded = math.ceil(max_barcodes / 100000) * 100000
    res = {}
    for batch in [0, 10, 50, 100, 500, 1000, 5000, 10000, 15000, 20000, 25000, max_rounded]:
        sub_df = check_df.loc[barcodes_df['Barcodes'] <= batch]
        check_df = check_df.loc[barcodes_df['Barcodes'] > batch]
        res[batch] = len(sub_df)
    fig = go.Figure()
    fig.add_trace(go.Bar(y=list(res.values()), x=[f'≤{i}' for i in res.keys()], text=list(res.values()), marker_color='navy'))
    fig.update_layout(template='simple_white', title='Barcode distribution')
    fig.update_xaxes(title='number of barcodes (categorized)')
    fig.update_yaxes(title='number of taxa')
    file_1 = Path('{}/{}.report_1.pdf'.format(output_directories[3], project_name))
    fig.write_image(file_1)


    ## 2) ranking distribution
    ratings = reference_db_df['rating'].values.tolist()
    counts = {i:ratings.count(i) for i in range(-20,51)}
    fig = go.Figure()
    fig.add_trace(go.Bar(y=list(counts.values()), x=list(counts.keys()), marker_color='navy'))
    fig.update_layout(template='simple_white', title='Rating distribution')
    fig.update_xaxes(title='rating', range=(-20,50))
    fig.update_yaxes(title='number of barcodes')
    fig.add_vrect(x0=9.5, x1=24.5, line_width=0, fillcolor="Peru", opacity=0.3, layer='below')
    fig.add_vrect(x0=24.5, x1=39.5, line_width=0, fillcolor="Silver", opacity=0.3, layer='below')
    fig.add_vrect(x0=39.5, x1=50.5, line_width=0, fillcolor="Gold", opacity=0.3, layer='below')
    file_2 = Path('{}/{}.report_2.pdf'.format(output_directories[3], project_name))
    fig.write_image(file_2)

    ## 3) database completeness
    test_taxon = 'family_name'
    taxa = sorted(reference_db_df[test_taxon].drop_duplicates().values.tolist())
    y_values = [reference_db_df[test_taxon].values.tolist().count(i) for i in taxa]
    height = optimal_height(y_values)
    fig = go.Figure()
    fig.add_trace(go.Bar(x=y_values[::-1], y=taxa[::-1], text=y_values[::-1], textposition='outside', cliponaxis=False, marker_color='navy', orientation='h'))
    fig.update_layout(template='simple_white',
                      width=1000,
                      height=height,
                      title = 'Barcodes per family'
                      )
    fig.update_xaxes(title='number of reference sequences')
    fig.update_yaxes(dtick='linear', automargin=True)
    file_3 = Path('{}/{}.report_3.pdf'.format(output_directories[3], project_name))
    fig.write_image(file_3)

    ## 4) percentage of phylogenetic stage
    taxa = sorted(reference_db_df[test_taxon].drop_duplicates().values.tolist())
    phylo_states = ['monophyletic', 'monophyletic (singleton)', 'paraphyletic', 'paraphyletic (LRT failed)', 'Ambiguous identification', 'Error']
    colors = ['Lightgreen', 'Green', 'Red', 'Darkred', 'Black', 'Grey']
    phylo_dict = {i:[] for i in phylo_states}
    for taxon in taxa:
        sub_df = reference_db_df.loc[reference_db_df[test_taxon] == taxon]
        n = len(sub_df)
        res = [sub_df['phylogeny'].values.tolist().count(i) for i in phylo_states]
        res_relative = [round(i/n*100,3) for i in res]
        for proportion, state in zip(res_relative, phylo_states):
            phylo_dict[state] = phylo_dict[state] + [proportion]

    fig = go.Figure()
    c = 0
    for state, y_values in phylo_dict.items():
        fig.add_trace(go.Bar(x=y_values[::-1], y=taxa[::-1], name=state, marker_color=colors[c], orientation='h'))
        c += 1
    fig.update_layout(template='simple_white',
                      width=1000,
                      height=height,
                      barmode='stack',
                      title = 'Phylogenetic approach'
                      )
    fig.update_xaxes(title='proportion of reference sequences')
    fig.update_yaxes(dtick='linear', automargin=True)
    file_4 = Path('{}/{}.report_4.pdf'.format(output_directories[3], project_name))
    fig.write_image(file_4)

    ## 5) sequence quality
    taxa = sorted(reference_db_df[test_taxon].drop_duplicates().values.tolist())
    quality_categories = ['Good', 'Sufficient', 'Bad']
    colors = ['Lightgreen', 'Green', 'Red']
    quality_dict = {i:[] for i in quality_categories}
    for taxon in taxa:
        sub_df = reference_db_df.loc[reference_db_df[test_taxon] == taxon]
        n = len(sub_df)
        res = [sub_df['sequence_quality'].values.tolist().count(i) for i in quality_categories]
        res_relative = [round(i/n*100,3) for i in res]
        for proportion, state in zip(res_relative, quality_dict):
            quality_dict[state] = quality_dict[state] + [proportion]

    fig = go.Figure()
    c = 0
    for state, y_values in quality_dict.items():
        fig.add_trace(go.Bar(x=y_values[::-1], y=taxa[::-1], name=state, marker_color=colors[c], orientation='h'))
        c += 1
    fig.update_layout(template='simple_white',
                      width=1000,
                      height=height,
                      barmode='stack',
                      title='Sequence quality'
                      )
    fig.update_xaxes(title='proportion of reference sequences')
    fig.update_yaxes(dtick='linear', automargin=True)
    file_5 = Path('{}/{}.report_5.pdf'.format(output_directories[3], project_name))
    fig.write_image(file_5)

    ## merge pdf files
    mergeFile = PyPDF2.PdfMerger()
    mergeFile.append(PyPDF2.PdfReader(file_1, 'rb'))
    mergeFile.append(PyPDF2.PdfReader(file_2, 'rb'))
    mergeFile.append(PyPDF2.PdfReader(file_3, 'rb'))
    mergeFile.append(PyPDF2.PdfReader(file_4, 'rb'))
    mergeFile.append(PyPDF2.PdfReader(file_5, 'rb'))
    merged_report = Path('{}/{}.report.pdf'.format(output_directories[3], project_name))
    mergeFile.write(merged_report)

    ## statistics table
    ## number of identifiers
    identifiers = reference_db_df['identification_by'].drop_duplicates()
    n_identifiers = len(identifiers)

    ## mean sequences per identifier
    sequence_per_identifier_dict = {i: len(reference_db_df.loc[reference_db_df['identification_by'] == i]) for i in identifiers}
    mean_sequences_per_identifier = np.mean(list(sequence_per_identifier_dict.values()))
    median_sequences_per_identifier = np.median(list(sequence_per_identifier_dict.values()))
    max_sequences_per_identifier = max(list(sequence_per_identifier_dict.values()))
    min_sequences_per_identifier = min(list(sequence_per_identifier_dict.values()))
    identifier_ranking = {}
    for i in [1,10,25,50,100,250,500,1000]:
        identifier_ranking[i] = len([values for key,values in sequence_per_identifier_dict.items() if values >= i])
    unique_species = len(reference_db_df['species_name'].drop_duplicates())
    phylogeny_counts = {}
    for i in sorted(reference_db_df['phylogeny'] .drop_duplicates().values.tolist()):
        phylogeny_counts[i] = len(reference_db_df.loc[reference_db_df['phylogeny'] == i])

########################################################################################################################

# settings_xlsx = '/Volumes/Coruscant/dbDNA/settings_mzb_mac.xlsx'

## load settings file
## collect user input from command line
if len(sys.argv) > 1:
    settings_xlsx = Path(sys.argv[1])

## otherwise set to default location
else:
    settings_xlsx = Path('./settings.xlsx')

## check if settings file is existing
if not os.path.isfile(settings_xlsx):
    user_input = input("Please provide the (full) PATH to a settings file:\n")
    settings_xlsx = Path(user_input)

## check if settings file is existing
if not os.path.isfile(settings_xlsx):
    print('Could not find the settings.xlsx!\n')
    print(settings_xlsx)

## run main script
else:
    # Collect tasks to run from settings file
    tasks = pd.read_excel(settings_xlsx, sheet_name='Tasks')
    data_source = tasks.values.tolist()[0][2]
    run_download = tasks.values.tolist()[1][2]
    run_extraction = tasks.values.tolist()[2][2]
    run_phylogeny = tasks.values.tolist()[3][2]
    run_rating = tasks.values.tolist()[4][2]
    run_create_database = tasks.values.tolist()[5][2]

    # Collect variables from settings file
    variables = pd.read_excel(settings_xlsx, sheet_name='Variables')
    project_name = variables.values.tolist()[0][2]
    taxa_list = variables.values.tolist()[1][2]
    identifier_whitelist = variables.values.tolist()[2][2]
    location_whitelist = variables.values.tolist()[3][2]
    output_folder = variables.values.tolist()[4][2]
    marker = variables.values.tolist()[5][2]
    min_rating = variables.values.tolist()[6][2]
    download_overwrite = variables.values.tolist()[7][2]
    alignment_overwrite = variables.values.tolist()[8][2]
    tree_overwrite = variables.values.tolist()[9][2]
    mafft_executable = variables.values.tolist()[10][2]
    iqtree_executable = variables.values.tolist()[11][2]
    mptp_executable = variables.values.tolist()[12][2]
    makeblastdb_exe = variables.values.tolist()[13][2]
    midori2_fasta = variables.values.tolist()[14][2]
    outgroup_fasta = variables.values.tolist()[15][2]
    lat_db = variables.values.tolist()[16][2]
    lon_db = variables.values.tolist()[17][2]
    d1 = variables.values.tolist()[18][2]
    d2 = variables.values.tolist()[19][2]
    d3 = variables.values.tolist()[20][2]
    use_coordinates = variables.values.tolist()[21][2]
    use_country = variables.values.tolist()[22][2]
    cpu_count = variables.values.tolist()[23][2]
    cpu_mode = variables.values.tolist()[24][2]

    ########################################################################################################################

    ## create output folders
    output_directories = create_new_project(project_name, output_folder)

    ## open the log file in append mode
    log_file = Path(f"{output_folder}/{project_name}.log")
    log_file = open(log_file, "a")

    ## create a custom stream that duplicates output to both console and log file
    class TeeStream:
        def __init__(self, *streams):
            self.streams = streams

        def write(self, data):
            for stream in self.streams:
                stream.write(data)

        def flush(self):
            for stream in self.streams:
                stream.flush()

    ## redirect stdout to both console and the log file
    sys.stdout = TeeStream(sys.stdout, log_file)

    ## test if enough cores are available
    available_cores = multiprocessing.cpu_count()
    if cpu_count > available_cores:
        cpu_count = multiprocessing.cpu_count() - 1
        print('{} - Not enough CPUs available. Defaulting to {} CPUs instead.'.format(datetime.now().strftime("%H:%M:%S"), cpu_count))
    else:
        print('{} - Using {} CPUs.'.format(datetime.now().strftime("%H:%M:%S"), cpu_count))

    ## print CPU usage
    if cpu_mode == 'Parallel':
        print('{} - {} mode: Using up to {} CPUs for running individual commands in parallel (recommended for smaller families).'.format(datetime.now().strftime("%H:%M:%S"), cpu_mode, cpu_count))
    elif cpu_mode == 'Multithreading':
        print('{} - {} mode: Using up to {} CPUs per individual command (better for larger families).'.format( datetime.now().strftime("%H:%M:%S"), cpu_mode, cpu_count))

    ## run scripts
    if run_download == 'yes':
        # BOLD systems workflow
        if data_source == 'BOLD':
            download_data_from_bold(taxa_list, output_directories, marker)

        # GenBank workflow
        elif data_source == 'NCBI':
            extract_MIDORI2_file(midori2_fasta, output_directories, taxa_list)

    if run_extraction == 'yes':
        # BOLD systems workflow
        if data_source == 'BOLD':
            extract_bold_json(output_directories, marker)

        # GenBank workflow
        elif data_source == 'NCBI':
            extract_genbank_files(output_directories)

    if run_phylogeny == 'yes':
        phylogenetic_approach(output_directories, mafft_executable, iqtree_executable, mptp_executable, cpu_count, cpu_mode)

    if run_rating == 'yes':
        rating_system(output_directories, identifier_whitelist, location_whitelist, project_name, cpu_count)

    if run_create_database == 'yes':
        create_database(output_directories, project_name, makeblastdb_exe)


    ## close the log file
    print('{} - Writing to log file...'.format(datetime.now().strftime("%H:%M:%S")))

    ## finish script
    print('\n{} - Done. Have a nice day!\n'.format(datetime.now().strftime("%H:%M:%S")))

