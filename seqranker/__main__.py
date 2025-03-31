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
import pickle
import multiprocessing
from itertools import chain
import subprocess
from datetime import datetime, timedelta
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
from collections import defaultdict
from io import StringIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Align.AlignInfo import SummaryInfo
from Bio import AlignIO
from Bio.Align import Alignment
from Bio.motifs import Motif
import os

####### STORE RUNTIMES #######

t0 = datetime.now()
print_lock = threading.Lock()

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

    dir_4 = Path('{}/{}'.format(output_directory, '5_report'))
    if not os.path.exists(dir_4):
        os.makedirs(dir_4)

    return [output_directory, dir_1, dir_2, dir_3, dir_4]

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
        family_table = family_dir / f"2_{family}_raw_barcodes.parquet.snappy"
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

####### RECORD BLACK LIST #######
def blacklist_filter(output_directories, record_blacklist):

    print(f'{datetime.now().strftime("%H:%M:%S")} - Removing records based on blacklist.')

    ## verify file
    if not os.path.isfile(record_blacklist):
        print(f'{datetime.now().strftime("%H:%M:%S")} - Could not find the provided blacklist file!')
        return

    ## read blacklist
    blacklist_df = pd.read_excel(record_blacklist).fillna('')
    families = set(blacklist_df['family_name'].values.tolist())
    print(f'{datetime.now().strftime("%H:%M:%S")} - Found {len(blacklist_df)} records from {len(families)} families to remove.')

    ## exclude records per family
    for family in families:
        family_table = Path(output_directories[2].joinpath(family, f'2_{family}_raw_barcodes.parquet.snappy'))
        records_to_remove = blacklist_df.loc[blacklist_df['family_name'] == family]['record_id'].values.tolist()
        if os.path.isfile(family_table):
            df = pd.read_parquet(family_table).fillna('')
            df_filtered = df[~df['record_id'].isin(records_to_remove)]
            df_filtered.to_parquet(family_table, compression='snappy')
            print(f'{datetime.now().strftime("%H:%M:%S")} - Removed {len(df) - len(df_filtered)} record(s) from {family}.')

    print(f'{datetime.now().strftime("%H:%M:%S")} - Finished removing records based on blacklist.\n')

####### PHYLOGENY #######

def dereplicate_sequences(fasta_file, fasta_file_derep, fasta_file_pkl, family):
    # Use defaultdict for simpler and faster appending
    res = defaultdict(list)
    for record in SeqIO.parse(fasta_file, "fasta"):
        res[str(record.seq)].append(record.id)

    # Calculate statistics
    total_reads = sum(len(ids) for ids in res.values())  # Total input reads
    total_dereplicated_reads = len(res)  # Total unique sequences
    dereplicated_percentage = round((total_dereplicated_reads / total_reads) * 100, 1)

    # Write dereplicated sequences to a new FASTA file
    i = 1
    dereplication_dict = {}
    with open(fasta_file_derep, 'w') as f:
        for seq, ids in res.items():
            header = f'ID_{i}'
            f.write(f'>{header}\n{seq}\n')
            dereplication_dict[header] = ids
            i += 1

    # Save dictionary to a file
    with open(fasta_file_pkl, 'wb') as f:
        pickle.dump(dereplication_dict, f)

    # Print statistics
    print(f'{datetime.now().strftime("%H:%M:%S")} - Dereplicated {family} from {total_reads} to {total_dereplicated_reads} reads ({dereplicated_percentage}%)')

## remove all empty files
def remove_empty_files(family, folder):

    family_dir = Path(f'{folder}/{family}')
    if os.path.isdir(family_dir):
        print(f'{datetime.now().strftime("%H:%M:%S")} - Cleaning up {family}.')
        for root, _, files in os.walk(family_dir):
            for file in files:
                file_path = Path(os.path.join(root, file))
                if os.path.isfile(file_path) and os.path.getsize(file_path) == 0:
                    os.remove(file_path)
                    print(f'  {datetime.now().strftime("%H:%M:%S")} - Removed {file_path.name}')

# Make sure only ACGTN are in the sequences
def convert_aln_file(aln_file, aln_file_adjusted):
    f = open(aln_file_adjusted,"w")
    for record in SeqIO.parse(aln_file, "fasta"):
        id = f">{record.id}\n"
        seq = ''.join([i if i in ['A', 'C', 'G', 'T', 'N'] else 'N' for i in str(record.seq).upper().replace('-', 'N')]) + '\n'
        if len(seq) != 0:
            f.write(id)
            f.write(seq)
    f.close()

def run_mafft(family, mafft_executable, folder, cpu_count):
    """
    Function to calculate multiple sequence alignment using MAFFT.
    """
    print(f'{datetime.now().strftime("%H:%M:%S")} - Starting MAFFT alignment for {family}.')

    # Subdirectories for output and temporary files
    family_dir = Path(f'{folder}/{family}')
    if not os.path.exists(family_dir):
        os.makedirs(family_dir)

    # family = 'Salifidae'

    # File paths
    fasta_file = Path(f'{family_dir}/1_{family}.fasta')
    fasta_file_derep = Path(f'{family_dir}/1_{family}_derep.fasta')
    fasta_file_pkl = Path(f'{family_dir}/1_{family}_ids.pkl')
    aln_file = Path(f'{family_dir}/1_{family}_alignment.fasta')
    aln_file_adjusted = Path(f'{family_dir}/1_{family}_alignment_adjusted.fasta')

    species_file_txt = Path(f'{family_dir}/3_{family}_species.txt')
    species_file_txt_snappy = Path(f'{family_dir}/3_{family}_species.parquet.snappy')
    family_table = Path(f'{family_dir}/2_{family}_raw_barcodes.parquet.snappy')
    ambiguous_barcodes_table = Path(f'{family_dir}/2_{family}_ambiguous_taxa.parquet.snappy')
    usable_barcodes_table = Path(f'{family_dir}/2_{family}_usable_taxa.parquet.snappy')
    df = pd.read_parquet(family_table)
    raw_records = df.values.tolist()

    print(f'{datetime.now().strftime("%H:%M:%S")} - Processing species records...')

    if not os.path.isfile(ambiguous_barcodes_table) and not os.path.isfile(usable_barcodes_table):
        # Extract records and categorize them
        record_ids = df['record_id'].values.tolist()
        usable_records = []
        ambiguous_assignments = []

        for record in tqdm(record_ids, desc=f'{datetime.now().strftime("%H:%M:%S")} - Collecting ambiguous and usable records'):
            sub_df = df.loc[df['record_id'] == record]
            species = sub_df['species'].values.tolist()[0]
            sequence = sub_df['nuc'].values.tolist()[0].upper().replace('-', '').strip('N')
            special_characters = '!@#$%^&*()-+?_=,.<>/\'\"0123456789'

            if any(c in special_characters for c in species):
                ambiguous_assignments.append([record, species, '', 'Ambiguous identification', '', -20]) # already write them in the final format for the rating
            else:
                usable_records.append([record, species, '', '', '', 0, sequence])

        # Save ambiguous taxa to file
        print(f'{datetime.now().strftime("%H:%M:%S")} - Writing ambiguous taxa to file...')
        ambiguous_table_df = pd.DataFrame(ambiguous_assignments, columns=['Sequence ID', 'Species Name', 'Species No.', 'Phylogeny', 'Clade', 'Rating'])
        ambiguous_table_df.to_parquet(ambiguous_barcodes_table, compression='snappy')

        # Save usable taxa to file
        print(f'{datetime.now().strftime("%H:%M:%S")} - Writing usable taxa to file...')
        usable_table_df = pd.DataFrame(usable_records, columns=['Sequence ID', 'Species Name', 'Cluster', 'Clade', 'State', 'Rating', 'Sequence'])
        usable_table_df.to_parquet(usable_barcodes_table, compression='snappy')
    else:
        print(f'{datetime.now().strftime("%H:%M:%S")} - Ambiguous and usable records were already stored.')
        print(f'{datetime.now().strftime("%H:%M:%S")} - Loading usable records...')
        usable_records = pd.read_parquet(usable_barcodes_table).values.tolist()

    # Determine whether alignment needs to be performed
    if len(usable_records) < 4:
        print(f'{datetime.now().strftime("%H:%M:%S")} - Skipped alignment for {family} (less than 4 sequences).')

        print(f'{datetime.now().strftime("%H:%M:%S")} - Creating dummy file for all records...')
        with open(species_file_txt, 'w') as f:
            for record in raw_records:
                f.write(f'{record[4]}__{record[21]}\n')

    # Check if analysis was already performed
    elif not os.path.isfile(aln_file) or not os.path.isfile(fasta_file) or not os.path.isfile(fasta_file_derep) or not os.path.isfile(fasta_file_pkl):

        print(f'{datetime.now().strftime("%H:%M:%S")} - Creating FASTA file...')
        with open(fasta_file, 'w') as f:
            for i, record in enumerate(usable_records):
                header = f'>{record[0]}__{record[1].replace(" ", "_")}__{i}\n'
                sequence = f'{record[-1]}\n'.upper().replace('-', '').strip('N')
                f.write(header + sequence)

        print(f'{datetime.now().strftime("%H:%M:%S")} - Dereplicating sequences...')
        dereplicate_sequences(fasta_file, fasta_file_derep, fasta_file_pkl, family)

        print(f'{datetime.now().strftime("%H:%M:%S")} - Running MAFFT alignment...')
        command = f"{mafft_executable} --auto --quiet --thread {cpu_count} --preservecase {fasta_file_derep} > {aln_file}"

        process = subprocess.Popen(command, shell=True)
        process.wait()

        print(f'{datetime.now().strftime("%H:%M:%S")} - Finished MAFFT alignment for {family}.')

    else:
        print(f'{datetime.now().strftime("%H:%M:%S")} - MAFFT alignment was already calculated for {family}.')

    print('')

## run species delimitation
def run_vsearch(family, vsearch_executable, folder, cpu_count):
    """
    Perform species delimitation for a given family using VSEARCH.
    """
    # family = 'Chironomidae'

    print(f'{datetime.now().strftime("%H:%M:%S")} - Starting species delimitation for {family}.')

    # Directory setup
    family_dir = Path(f'{folder}/{family}')
    os.makedirs(family_dir, exist_ok=True)

    clusters_folder = Path(f'{family_dir}/clusters/')
    os.makedirs(clusters_folder, exist_ok=True)

    # File paths
    fasta_file_pkl = Path(f'{family_dir}/1_{family}_ids.pkl')
    aln_file = Path(f'{family_dir}/1_{family}_alignment.fasta')
    aln_file_adjusted = Path(f'{family_dir}/1_{family}_alignment_adjusted.fasta')
    centroids_file = Path(f'{family_dir}/1_{family}.centroids')
    species_file_txt_snappy = Path(f'{family_dir}/3_{family}_species.parquet.snappy')
    family_table = Path(f'{family_dir}/2_{family}_raw_barcodes.parquet.snappy')

    # Load raw records
    df = pd.read_parquet(family_table)
    raw_records = df.values.tolist()

    # Check if analysis is already done or alignment is missing
    if os.path.isfile(species_file_txt_snappy):
        print(f'{datetime.now().strftime("%H:%M:%S")} - {family} was already analysed!')
        print('')
        return

    if not os.path.isfile(aln_file) and not os.path.isfile(aln_file_adjusted):
        print(f'{datetime.now().strftime("%H:%M:%S")} - {family} is missing an alignment file. Cannot calculate clusters!')

        species_data = [[record[4], record[21], '', 'Error: Vsearch could not find alignment file', '', 0] for record in raw_records]
        species_df = pd.DataFrame(
            species_data, columns=['Sequence ID', 'Species Name', 'Species No.', 'Phylogeny', 'Clade', 'Rating']
        )
        species_df = species_df.astype('string')
        species_df.to_parquet(species_file_txt_snappy, index=False, compression='snappy')
        print('')
        return

    print(f'{datetime.now().strftime("%H:%M:%S")} - Preparing alignment for VSEARCH.')
    convert_aln_file(aln_file, aln_file_adjusted)
    print(f'{datetime.now().strftime("%H:%M:%S")} - Adjusted alignment for {family}.')

    # Run VSEARCH for species delimitation
    print(f'{datetime.now().strftime("%H:%M:%S")} - Running VSEARCH for {family}.')

    # Collect clusters
    clusters = glob.glob(str(clusters_folder.joinpath(f'{family}_*')))
    if not os.path.isfile(centroids_file) or not clusters:

        clusters_path = Path(f'{clusters_folder}/{family}_')
        # Clean up previous cluster files
        for file in glob.glob(f"{clusters_folder}/*"):
            os.remove(file)

        command = (
            f"{vsearch_executable} --cluster_size {aln_file_adjusted} "
            f"--id 0.99 --threads {cpu_count} --centroids {centroids_file} --clusters {str(clusters_path)}"
        )
        subprocess.run(command, shell=True)

    else:
        print(f'{datetime.now().strftime("%H:%M:%S")} - {family} Clusters were already calculated.')

    # Update clusters
    clusters = glob.glob(str(clusters_folder.joinpath(f'{family}_*')))
    if not clusters:
        print(f'{datetime.now().strftime("%H:%M:%S")} - No clusters found for {family}.')
        species_data = [[record[4], record[21], '', 'Error: Vsearch was unable to produce clusters', '', 0] for record in raw_records]
        species_df = pd.DataFrame(
            species_data, columns=['Sequence ID', 'Species Name', 'Species No.', 'Phylogeny', 'Clade', 'Rating']
        )
        species_df = species_df.astype('string')
        species_df.to_parquet(species_file_txt_snappy, index=False, compression='snappy')
        print(f'{datetime.now().strftime("%H:%M:%S")} - Finished {family}, but failed species delimitation.')
        return

    # Load dereplication dictionary
    with open(fasta_file_pkl, 'rb') as f:
        dereplication_dict = pickle.load(f)

    res = []
    added_ids = []  # Track added IDs
    for i, file in enumerate(clusters):
        cluster = [record.id for record in SeqIO.parse(file, "fasta")]
        added_ids.extend(cluster)
        ids = [dereplication_dict[id] for id in cluster]
        flattened_ids = list(chain.from_iterable(ids))
        rows = [id.split('__') for id in flattened_ids]
        cluster_df = pd.DataFrame(rows, columns=['Sequence ID', 'Species Name', 'No.'])

        n_records = len(cluster_df)
        n_species = len(set(cluster_df['Species Name']))

        clade = f'{family}_{i + 1}'
        if n_records >= 2 and n_species == 1:
            phylogeny = 'monophyletic'
        elif n_records < 2 and n_species == 1:
            phylogeny = 'monophyletic (singleton)'
        elif n_species != 1:
            phylogeny = 'polyphyletic'
        else:
            phylogeny = 'polyphyletic (insufficient data)'

        for record in cluster_df.values.tolist():
            res.append(record + [phylogeny, clade, 0])

    # Add IDs removed based on overlap matrix
    for id in dereplication_dict.keys():
        if id not in added_ids:
            ids = dereplication_dict[id]
            rows = [id.split('__') for id in ids]
            cluster_df = pd.DataFrame(rows, columns=['Sequence ID', 'Species Name', 'No.'])
            for record in cluster_df.values.tolist():
                res.append(record + ['spurious sequence', '', 0])

    # Write species delimitation table
    df_out = pd.DataFrame(res, columns=['Sequence ID', 'Species Name', 'Species No.', 'Phylogeny', 'Clade', 'Rating'])
    df_out = df_out.astype('string')  # Ensure compatibility for Parquet
    df_out.to_parquet(species_file_txt_snappy, index=False, compression='snappy')
    print(f'{datetime.now().strftime("%H:%M:%S")} - Finished species delimitation for {family}.')

    # Handle missing species file
    if not os.path.isfile(species_file_txt_snappy):
        species_data = [[record[4], record[21], '', 'Error', '', 0] for record in raw_records]
        species_df = pd.DataFrame(
            species_data, columns=['Sequence ID', 'Species Name', 'Species No.', 'Phylogeny', 'Clade', 'Rating']
        )
        species_df = species_df.astype('string')
        species_df.to_parquet(species_file_txt_snappy, index=False, compression='snappy')
        print(f'{datetime.now().strftime("%H:%M:%S")} - Finished {family}, but failed species delimitation.')

    print('')

# main script for phylogenetic approach
def phylogenetic_approach(output_directories, mafft_executable, vsearch_executable, similarity_threshold, cpu_count):

    print('{} - Starting phylogenetic approach.'.format(datetime.now().strftime("%H:%M:%S")))

    # Store fasta, alignments and tree
    folder = output_directories[2]

    print('{} - Sorting families by size and starting with the smallest.'.format(datetime.now().strftime("%H:%M:%S")))

    # Extract families
    files = glob.glob('{}/*/*_raw_barcodes.parquet.snappy'.format(output_directories[2]))
    # Sort them by size
    files = {file: size for file, size in sorted({file: os.path.getsize(file) for file in files}.items(), key=lambda item: item[1], reverse=False)}
    # Collect family names
    families = [Path(i).name.replace('_raw_barcodes.parquet.snappy', '').replace('2_', '') for i in files.keys()]

    for family, size in files.items():
        fam = Path(family).name.replace('_raw_barcodes.parquet.snappy', '').replace('2_', '')
        size_mb = round(size / (1024 ** 2), 2)
        print(f'{fam} ({size_mb} MB)')
    print('')

    print('{} - Using {}/{} CPUs in Multithreading mode.'.format(datetime.now().strftime("%H:%M:%S"), cpu_count, cpu_count+1))
    print('{} - Starting species delimitation.'.format(datetime.now().strftime("%H:%M:%S")))

    print('{} - Starting pre-cleanup to remove empty files.'.format(datetime.now().strftime("%H:%M:%S")))
    [remove_empty_files(family, folder) for family in families]
    print('{} - Finished pre-cleanup.\n'.format(datetime.now().strftime("%H:%M:%S")))

    ## Mafft alignment
    print('{} - Starting mafft alignments.'.format(datetime.now().strftime("%H:%M:%S")))
    [run_mafft(family, mafft_executable, folder, cpu_count) for family in families]
    print('{} - Finished mafft alignments.\n'.format(datetime.now().strftime("%H:%M:%S")))

    # Vsearch clustering
    print('{} - Starting vsearch clustering.'.format(datetime.now().strftime("%H:%M:%S")))
    [run_vsearch(family, vsearch_executable, folder, cpu_count) for family in families]
    print('{} - Finished vsearch clustering.\n'.format(datetime.now().strftime("%H:%M:%S")))

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

####### RATING #######

# function to rate each family
def rate_family(output_directories, identifier_whitelist_lst, location_whitelist_lst, family, use_coordinates, use_country, d1, d2, d3):

    # family = 'Arthropleidae'
    # check_coordinates = True
    # df  = pd.read_parquet('/Volumes/Coruscant/dbDNA/Small_db_BarCodeBank/2_phylogeny/Viviparidae/Viviparidae_raw_barcodes.parquet.snappy')
    # voucher_status
    # identification_method

    # Sub directories for output and temporary files
    family_dir = Path('{}/{}'.format(output_directories[2], family))
    if not os.path.exists(family_dir):
        os.makedirs(family_dir)

    ## relevant files
    species_file_txt_snappy = Path(f'{family_dir}/3_{family}_species.parquet.snappy')
    ambiguous_taxa_table = Path(f'{family_dir}/2_{family}_ambiguous_taxa.parquet.snappy')
    ratings_file = Path(f'{family_dir}/4_{family}_ratings.parquet.snappy')
    family_table = Path(f'{family_dir}/2_{family}_raw_barcodes.parquet.snappy')

    if os.path.isfile(ratings_file):
        print('{} - Rating was already conducted for {}.'.format(datetime.now().strftime("%H:%M:%S"), family))
        return

    # Collect information about records
    df1 = pd.read_parquet(family_table).fillna('')

    ## store ratings in list
    all_ratings_list = []

    # Location white list
    main_country = [i for i in location_whitelist_lst['main country'].values.tolist() if i != '']
    neighbour_countries = [i for i in location_whitelist_lst['neighbour_countries'].values.tolist() if i != '']
    continent = [i for i in location_whitelist_lst['continent'].values.tolist() if i != '']

    ## read result files
    df2 = pd.read_parquet(species_file_txt_snappy).fillna('')
    df3 = pd.read_parquet(ambiguous_taxa_table).fillna('')

    # Remove empty or all-NaN columns from df2 and df3
    df2 = df2.dropna(axis=1, how='all')
    df3 = df3.dropna(axis=1, how='all')

    # Merge clustered and excluded (ambiguous) data
    merged_df = pd.concat([df2, df3], ignore_index=True)

    ## rate each individual barcode of the family
    for record in merged_df.values.tolist():
        ## collect information about the record
        # family file
        record_id = record[0]
        species_group = record[4]
        phylogeny = record[3]
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
        image_urls = '' # not available from boldsystems v5
        markercode = raw_record['marker_code'].values.tolist()[0]
        nucleotides = raw_record['nuc'].values.tolist()[0].upper().replace('-', '').strip('N')
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
            rating += 6
            sequence_quality = 'Good'
            rating_calc.append('+6 (sequence quality)')

        ## bad sequence quality
        not_allowed_chars = [i for i in set(nucleotides_trimmed) if i not in allowed_chars]
        if len(not_allowed_chars) != 0:
            n_not_allowed_chars = sum([nucleotides_trimmed.count(i) for i in not_allowed_chars])
            rel = n_not_allowed_chars / len(nucleotides_trimmed) * 100
            if rel >= 2:
                rating -= 10
                sequence_quality = 'Bad'
                rating_calc.append('-10 (sequence quality)')

        ## sequence length (>= 500 bp are accepted as barcode)
        if len(nucleotides_trimmed) >= 500:
            rating += 5
            rating_calc.append('+5 (barcode length)')

        ## reverse BIN taxonomy
        keywords = ["bin", "bold", "tree"]
        exclude = ["morphology", "bold:"]
        value = identification_method.lower()
        if any(word in value for word in keywords) and not any(word in value for word in exclude):
            rating -= 10
            reverse_bin = 'Investigate'
            rating_calc.append('-10 (rev. BIN identification)')

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
    ratings_df.to_parquet(ratings_file, index=False, compression='snappy')

    print('{} - Finished rating for {}.'.format(datetime.now().strftime("%H:%M:%S"), family))

# main script for rating algorithm
def rating_system(output_directories, identifier_whitelist, location_whitelist, project_name, cpu_count):

    print('{} - Collecting raw records.'.format(datetime.now().strftime("%H:%M:%S")))

    # Extract families
    files = glob.glob('{}/*/*_raw_barcodes.parquet.snappy'.format(output_directories[2]))
    families = sorted([Path(i).name.replace('_raw_barcodes.parquet.snappy', '').replace('2_', '') for i in files])

    identifier_whitelist_lst = sorted(set(pd.read_excel(identifier_whitelist, sheet_name='Identifier_Whitelist').fillna('')['Name (as used on BOLD)'].values.tolist()))
    location_whitelist_lst = pd.read_excel(location_whitelist).fillna('')

    ## for testing
    # files_to_delete = glob.glob('/Volumes/Coruscant/dbDNA/FEI_genera_BarCodeBank/2_phylogeny/*/*.ratings.parquet.snappy')
    # [os.remove(file) for file in files_to_delete if os.path.isfile(file)]

    print('{} - Starting to rate sequences.'.format(datetime.now().strftime("%H:%M:%S")))

    ## rate all families in parallel
    Parallel(n_jobs=cpu_count, backend='threading')(delayed(rate_family)(output_directories, identifier_whitelist_lst, location_whitelist_lst, family, use_coordinates, use_country, d1, d2, d3) for family in families)

    print('{} - Finished rating for all families.'.format(datetime.now().strftime("%H:%M:%S")))

    # Initialize an empty list to hold all dataframes
    dfs = []

    print('{} - Collecting ratings for all families.'.format(datetime.now().strftime("%H:%M:%S")))
    for family in tqdm(families):
        family_dir = Path('{}/{}'.format(output_directories[2], family))
        ratings_snappy = Path(f'{family_dir}/4_{family}_ratings.parquet.snappy')

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
    ratings_df.to_parquet(output_file_1, index=False, compression='snappy')
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

####### BLAST DATABASE #######

def create_database(output_directories, project_name, makeblastdb_exe):

    print('{}: Starting to build a new database.'.format(datetime.now().strftime('%H:%M:%S')))

    database_snappy = Path(f'{output_directories[3]}/{project_name}.BarCodeBank.parquet.snappy')
    df = pd.read_parquet(database_snappy)

    # Create fasta file
    fasta_file = Path(f'{output_directories[3]}/{project_name}.BarCodeBank.fasta')

    records = []
    for sequence in df.values.tolist():
        nucleotides = sequence[-2].upper().replace('-', '').strip('N')
        seq_id = sequence[1]
        taxonomy = '__'.join([i.replace(" ", "_") for i in sequence[3:9]])
        header = f'{seq_id}__{taxonomy}'
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

####### REPORT #######

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

def get_higher_taxon_from_gbif(family_name, higher_taxon):
    # GBIF API URL for taxon search
    url = f"https://api.gbif.org/v1/species/match?name={family_name}&rank=family"

    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()
        if higher_taxon in data and higher_taxon + 'Key' in data:
            return data[higher_taxon], data[higher_taxon + 'Key']
        else:
            return "f{higher_taxon} not found", None
    else:
        return f"Error: {response.status_code}", None

def create_report(output_directories, project_name, taxa_list):

    print(f'\n{datetime.now().strftime("%H:%M:%S")} - Creating report for project: {project_name}.')

    ## import files
    barcode_bank_file = Path('{}/{}.BarCodeBank.parquet.snappy'.format(output_directories[3], project_name))
    reference_db_df = pd.read_parquet(barcode_bank_file).fillna('')

    # Extract families
    dfs = []
    files = glob.glob('{}/*/*_raw_barcodes.parquet.snappy'.format(output_directories[2]))
    families = sorted([Path(i).name.replace('_raw_barcodes.parquet.snappy', '').replace('2_', '') for i in files])
    for family in families:
        family_dir = Path('{}/{}'.format(output_directories[2], family))
        ratings_snappy = family_dir / '2_{}_raw_barcodes.parquet.snappy'.format(family)
        if not ratings_snappy.is_file():
            print('{} - WARNING: No ratings file found for {}!!'.format(datetime.now().strftime("%H:%M:%S"), family))
        else:
            df = pd.read_parquet(ratings_snappy)
            dfs.append(df)  # Append the dataframe to the list
    # Concatenate all dataframes in the list
    raw_data_df = pd.concat(dfs, ignore_index=True)

    ## check for missing reference barcodes
    ids_0 = raw_data_df['record_id'].drop_duplicates()
    ids_1 = reference_db_df['record_id'].drop_duplicates()
    len_0 = len(ids_0)
    len_1 = len(ids_1)
    n_shared = len(set(ids_0) & set(ids_1))
    shared_perc = n_shared / len_0 * 100
    print(f'{datetime.now().strftime("%H:%M:%S")} - Detected reference barcodes: {shared_perc} %, {n_shared} in total.')

    ## 1) input taxa without barcodes
    print(f'{datetime.now().strftime("%H:%M:%S")} - Collecting number of records per input input taxon.')
    taxa_list_df = pd.read_excel(taxa_list)
    all_levels = reference_db_df[['order_name', 'class_name', 'family_name', 'genus_name', 'species_name']].values.tolist()
    flattened_levels = list([item for sublist in all_levels for item in sublist if item != ''])
    res = {}
    for taxon in taxa_list_df.values.tolist():
        taxon = taxon[0]
        res[taxon] = flattened_levels.count(taxon)
    barcodes_df = pd.DataFrame([[i,j] for i,j in res.items()], columns=['Taxon', 'Barcodes']).sort_values('Barcodes')
    barcodes_df.to_excel(Path('{}/{}.report_barcodes.xlsx'.format(output_directories[4], project_name)), index=False)

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
    file_1 = Path('{}/{}.report_1.pdf'.format(output_directories[4], project_name))
    fig.write_image(file_1)

    ## 2) ranking distribution
    print(f'{datetime.now().strftime("%H:%M:%S")} - Calculating ranking distribution.')
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
    file_2 = Path('{}/{}.report_2.pdf'.format(output_directories[4], project_name))
    fig.write_image(file_2)

    ## 3) database completeness
    print(f'{datetime.now().strftime("%H:%M:%S")} - Calculated database completeness.')
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
    file_3 = Path('{}/{}.report_3.pdf'.format(output_directories[4], project_name))
    fig.write_image(file_3)

    ## 4) percentage of phylogenetic stage
    print(f'{datetime.now().strftime("%H:%M:%S")} - Assessing phylogenetic states.')
    taxa = sorted(reference_db_df[test_taxon].drop_duplicates().values.tolist())
    phylo_states = ['monophyletic', 'monophyletic (singleton)', 'polyphyletic', 'Ambiguous identification', 'Error', 'spurious sequence']
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
    file_4 = Path('{}/{}.report_4.pdf'.format(output_directories[4], project_name))
    fig.write_image(file_4)

    ## 5) sequence quality
    print(f'{datetime.now().strftime("%H:%M:%S")} - Assessing sequence quality.')
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
    file_5 = Path('{}/{}.report_5.pdf'.format(output_directories[4], project_name))
    fig.write_image(file_5)

    ## 6) Calculate barcode coverage

    print(f'{datetime.now().strftime("%H:%M:%S")} - Assessing barcode coverage.')
    reference_xlsx = Path(reference_xlsx)
    if not reference_xlsx.exists():
        print(f'{datetime.now().strftime("%H:%M:%S")} - Could not find a reference list!')
    else:
        reference_df = pd.read_excel(reference_xlsx).fillna('')
        reference_species = []
        higher_taxon = reference_taxon

        if 'Species' in reference_df.columns.tolist():
            dropped_species = []
            for h_taxon, species in reference_df[[higher_taxon, 'Species']].values.tolist():
                if species != '' and all(char.isalpha() or char.isspace() for char in species):
                    name = f"{species.split(' ')[0]} {species.split(' ')[1]}"
                    reference_species.append([h_taxon.title(), name])
                else:
                    dropped_species.append(species)

            ref_taxa = pd.DataFrame(reference_species, columns=[higher_taxon, 'species_name'])
            db_species = sorted(reference_db_df['species_name'].drop_duplicates().values.tolist())
            db_genera = sorted(set([i.split(' ')[0] for i in db_species]))

            # print dropped species
            rel_dropped = round(len(dropped_species) / (len(dropped_species) + len(reference_species)) * 100, 2)
            print(f'{datetime.now().strftime("%H:%M:%S")} - Dropped {len(dropped_species)} ({rel_dropped}%) taxa that did not meet species recognition criteria!')

            barcode_coverage_dict = {}
            barcode_coverage_list = []
            for h_taxon in sorted(ref_taxa[higher_taxon].drop_duplicates().values.tolist(), reverse=True):
                h_taxon_species = ref_taxa[ref_taxa[higher_taxon] == h_taxon]['species_name'].drop_duplicates().values.tolist()
                species_present = []
                genus_present = []
                missing = []
                for species in h_taxon_species:
                    genus = species.split(' ')[0]
                    if species in db_species:
                        species_present.append(species)
                        res = [h_taxon, species, 'species present']
                    elif genus in db_genera:
                        genus_present.append(species)
                        res = [h_taxon, species, 'genus present']
                    else:
                        missing.append(species)
                        res = [h_taxon, species, 'species absent']
                    barcode_coverage_list.append(res)
                total = len(species_present) + len(genus_present) + len(missing)
                barcode_coverage_dict[h_taxon] = [len(species_present)/total, len(genus_present)/total, len(missing)/total, total]

            # Barcode Coverage
            fig = go.Figure()
            species_present = [i[0] for i in barcode_coverage_dict.values()]
            genus_present = [i[1] for i in barcode_coverage_dict.values()]
            total = [i[3] for i in barcode_coverage_dict.values()]
            families = [i for i in barcode_coverage_dict.keys()]
            fig.add_trace(go.Bar(x=species_present, y=families, name='species present', marker_color='Green', orientation='h'))
            fig.add_trace(go.Bar(x=genus_present, y=families, name='genus present', marker_color='Lightgreen', orientation='h'))
            # Add annotations for total counts
            for fam, count, species in zip(families, total, species_present):
                fig.add_annotation(
                    x=1.05,  # Position the text at the middle of the species bar
                    y=fam,
                    text=str(count),  # Convert total count to string
                    showarrow=False,  # No arrow, just text
                    font=dict(size=12, color='black'),
                    xanchor="center",  # Align text in the middle
                    yanchor="middle"
                )
            fig.update_layout(template='simple_white',
                              width=1000,
                              height= optimal_height(families),
                              barmode='stack',
                              title=f'Barcode coverage: "{higher_taxon}"'
                              )
            fig.update_xaxes(title='Barcode coverage (species)', dtick=0.1)
            fig.update_yaxes(dtick='linear', automargin=True)
            file_6 = Path('{}/{}.report_6.pdf'.format(output_directories[4], project_name))
            fig.write_image(file_6)

            # also write the table
            barcode_coverage_df = pd.DataFrame(barcode_coverage_list, columns=[higher_taxon, 'Species', 'Status'])
            barcode_coverage_xlsx = Path('{}/{}.barcode_coverage.xlsx'.format(output_directories[4], project_name))
            barcode_coverage_df.to_excel(barcode_coverage_xlsx, index=False)

        else:
            print(f'{datetime.now().strftime("%H:%M:%S")} - Could not find "{higher_taxon}" in the reference list!')

        ## 7) rating per taxonomic level
        print(f'{datetime.now().strftime("%H:%M:%S")} - Assessing rating per taxonomic level.')
        levels = ['phylum_name', 'class_name', 'order_name', 'family_name']
        i = 1
        for level in tqdm(levels):
            fig = go.Figure()
            taxa = sorted(reference_db_df[level].drop_duplicates().values.tolist(), reverse=True)
            res = {taxon:reference_db_df[reference_db_df[level]==taxon]['rating'].values.tolist() for taxon in taxa}
            for taxon in taxa:
                fig.add_trace(go.Box(x=res[taxon], y=[taxon]*len(res[taxon]), line_width=1, marker_size=1, marker_color='Navy', orientation='h'))
                fig.add_annotation(text=f'Ø{str(round(np.mean(res[taxon]), 2))}', x=51, y=taxon, xanchor="left", showarrow=False, font=dict(size=8, color='black'))
            fig.update_layout(template='simple_white',
                              width=1000,
                              height= optimal_height(taxa),
                              title=f'Rating distribution: "{level}"',
                              showlegend=False
                              )
            fig.update_xaxes(title=f'Rating distribution', range=(-21,51), dtick=10)
            fig.update_yaxes(dtick='linear')
            fig.add_vrect(x0=9.5, x1=24.5, line_width=0, fillcolor="Peru", opacity=0.3, layer='below')
            fig.add_vrect(x0=24.5, x1=39.5, line_width=0, fillcolor="Silver", opacity=0.3, layer='below')
            fig.add_vrect(x0=39.5, x1=50.5, line_width=0, fillcolor="Gold", opacity=0.3, layer='below')
            file = Path('{}/{}.report_{}_7.pdf'.format(output_directories[4], project_name, i))
            fig.write_image(file)
            i += 1

    ## merge pdf files
    print(f'{datetime.now().strftime("%H:%M:%S")} - Creating report file.')
    mergeFile = PyPDF2.PdfMerger()
    mergeFile.append(PyPDF2.PdfReader(file_1, 'rb'))
    os.remove(file_1)
    mergeFile.append(PyPDF2.PdfReader(file_2, 'rb'))
    os.remove(file_2)
    mergeFile.append(PyPDF2.PdfReader(file_3, 'rb'))
    os.remove(file_3)
    mergeFile.append(PyPDF2.PdfReader(file_4, 'rb'))
    os.remove(file_4)
    mergeFile.append(PyPDF2.PdfReader(file_5, 'rb'))
    os.remove(file_5)
    if file_6.is_file():
        mergeFile.append(PyPDF2.PdfReader(file_6, 'rb'))
        os.remove(file_6)
    file_7 = glob.glob(str(Path('{}/{}*_7.pdf'.format(output_directories[4], project_name))))
    for file in file_7:
        mergeFile.append(PyPDF2.PdfReader(file, 'rb'))
        os.remove(file)
    merged_report = Path('{}/{}.report.pdf'.format(output_directories[4], project_name))
    mergeFile.write(merged_report)

    print(f'{datetime.now().strftime("%H:%M:%S")} - Finished creating report for: {project_name}.\n')

def db_comparison(output_directories, project_name, version, project_name_or):
    if version == 1:
        print(f'{datetime.now().strftime("%H:%M:%S")} - Database comparison requires at least version 2.')
        return

    latest_db = Path(f"{output_directories[3]}/{project_name}.BarCodeBank.parquet.snappy")
    latest_df = pd.read_parquet(latest_db).fillna('')
    current_version = f'_v{version}'
    expected_versions = [i + 1 for i in range(0, version - 1)]
    old_dbs = [str(latest_db).replace(current_version, f'_v{i}') for i in expected_versions]

    for file in old_dbs:
        if not os.path.isfile(file):
            print(f'{datetime.now().strftime("%H:%M:%S")} - Missing file: {file}')
            return

    ratings_summary_list = []  # List to collect ratings summaries
    records_comparison_list = []  # List to collect records comparison results
    species_comparison_list = []  # List to collect species comparison results

    for test_version, file in enumerate(old_dbs + [latest_db], start=1):
        comparison_df = pd.read_parquet(file).fillna('')

        if file != str(latest_db):  # Skip unnecessary comparisons for the latest version
            # Records comparison
            new_set = set(latest_df['record_id'].values.tolist())
            old_set = set(comparison_df['record_id'].values.tolist())

            n_shared = len(new_set & old_set)
            n_old_only = len(old_set - new_set)
            n_new_only = len(new_set - old_set)

            records_comparison_list.append({
                'Version': test_version,
                'Shared Records': n_shared,
                'Exclusive to Version': n_old_only,
                'Exclusive to Latest': n_new_only
            })

            # Species comparison
            new_species_set = set(latest_df['species_name'].values.tolist())
            old_species_set = set(comparison_df['species_name'].values.tolist())

            n_shared_species = len(new_species_set & old_species_set)
            n_old_only_species = len(old_species_set - new_species_set)
            n_new_only_species = len(new_species_set - old_species_set)

            species_comparison_list.append({
                'Version': test_version,
                'Shared Species': n_shared_species,
                'Exclusive Species to Version': n_old_only_species,
                'Exclusive Species to Latest': n_new_only_species
            })

        # Ratings comparison (absolute numbers)
        gold = (comparison_df['rating'] >= 40).sum()
        silver = (comparison_df['rating'] >= 25).sum()
        bronze = (comparison_df['rating'] >= 10).sum()
        unreliable = (comparison_df['rating'] < 10).sum()

        total = gold + silver + bronze + unreliable

        # Calculate relative proportions
        gold_pct = round((gold / total) * 100, 2) if total > 0 else 0
        silver_pct = round((silver / total) * 100, 2) if total > 0 else 0
        bronze_pct = round((bronze / total) * 100, 2) if total > 0 else 0
        unreliable_pct = round((unreliable / total) * 100, 2) if total > 0 else 0

        # Collect both absolute and relative ratings summaries into the list
        ratings_summary_list.append({
            'Version': test_version if file != str(latest_db) else f'Latest (v{version})',
            'Gold': gold,
            'Silver': silver,
            'Bronze': bronze,
            'Unreliable': unreliable,
            'Gold (%)': gold_pct,
            'Silver (%)': silver_pct,
            'Bronze (%)': bronze_pct,
            'Unreliable (%)': unreliable_pct
        })

    # Convert lists to DataFrames
    # Ratings df
    ratings_summary_df = pd.DataFrame(ratings_summary_list)
    fig = go.Figure()
    ranks = ['Gold', 'Silver', 'Bronze', 'Unreliable'][::-1]
    colors = ['Gold', 'Silver', 'Peru', 'White'][::-1]
    for c, rank in enumerate(ranks):
        sub_df = ratings_summary_df[['Version', rank]]
        y_values = sub_df[rank].values.tolist()
        x_values = [f'v{i}' for i in sub_df['Version'].values.tolist()]
        fig.add_trace(go.Bar(y=y_values, x=x_values, name=rank, text=[str(i) if i != 0 else '' for i in y_values], textposition='outside', marker_line_color='black', marker_color=colors[c]))
    fig.update_layout(barmode='stack', template='simple_white')
    fig.update_yaxes(title='# records')
    file_1 = output_directories[4].joinpath('comparison_1.pdf')
    fig.write_image(file_1)

    records_comparison_df = pd.DataFrame(records_comparison_list)
    fig = go.Figure()
    ranks = ['Shared Records', 'Exclusive to Version', 'Exclusive to Latest']
    colors = ['Gold', 'Silver', 'Peru']
    for c, rank in enumerate(ranks):
        sub_df = records_comparison_df[['Version', rank]]
        y_values = sub_df[rank].values.tolist()
        x_values = [f'v{i}' for i in sub_df['Version'].values.tolist()]
        fig.add_trace(go.Bar(y=y_values, x=x_values, name=rank, text=[str(i) if i != 0 else '' for i in y_values], textposition='outside', marker_line_color='black', marker_color=colors[c]))
    fig.update_layout(barmode='stack', template='simple_white')
    fig.update_yaxes(title='# records')
    file_2 = output_directories[4].joinpath('comparison_2.pdf')
    fig.write_image(file_2)

    species_comparison_df = pd.DataFrame(species_comparison_list)
    fig = go.Figure()
    ranks = ['Shared Species', 'Exclusive Species to Version', 'Exclusive Species to Latest']
    colors = ['Gold', 'Silver', 'Peru']
    for c, rank in enumerate(ranks):
        sub_df = species_comparison_df[['Version', rank]]
        y_values = sub_df[rank].values.tolist()
        x_values = [f'v{i}' for i in sub_df['Version'].values.tolist()]
        fig.add_trace(go.Bar(y=y_values, x=x_values, name=rank, text=[str(i) if i != 0 else '' for i in y_values], textposition='outside', marker_line_color='black', marker_color=colors[c]))
    fig.update_layout(barmode='stack', template='simple_white')
    fig.update_yaxes(title='# records')
    file_3 = output_directories[4].joinpath('comparison_3.pdf')
    fig.write_image(file_3)

    ## merge pdf files
    print(f'{datetime.now().strftime("%H:%M:%S")} - Creating summary file.')
    mergeFile = PyPDF2.PdfMerger()
    mergeFile.append(PyPDF2.PdfReader(file_1, 'rb'))
    os.remove(file_1)
    mergeFile.append(PyPDF2.PdfReader(file_2, 'rb'))
    os.remove(file_2)
    mergeFile.append(PyPDF2.PdfReader(file_3, 'rb'))
    os.remove(file_3)
    merged_report = Path('{}/{}_version_comparison.pdf'.format(output_directories[4], project_name))
    mergeFile.write(merged_report)

    print(f'{datetime.now().strftime("%H:%M:%S")} - Finished version comparison: {project_name}.\n')

####### MISC #######

def filter_BTL_taxa():
    OTL_xlsx = '/Users/tillmacher/Desktop/TTT_projects/GeDNA_MZB_20_21/perlodes_TTT_conversion.xlsx'
    OTL_df = pd.read_excel(OTL_xlsx).fillna('')
    OTL_species = []
    for species in OTL_df['Species'].values.tolist():
        if species != '' and all(char.isalnum() or char.isspace() for char in species):
            name = f"{species.split(' ')[0]} {species.split(' ')[1]}"
            OTL_species.append(name)

    dbDNA_parquet = '/Volumes/Coruscant/dbDNA/MZB_Germany_v1_BarCodeBank/3_BarCodeBank/MZB_Germany_v1.BarCodeBank.parquet.snappy'
    dbDNA_df = pd.read_parquet(dbDNA_parquet)
    dbDNA_df_OTL_groups = sorted(set(dbDNA_df.loc[dbDNA_df['species_name'].isin(OTL_species)]['species_group'].values.tolist()))
    if '' in dbDNA_df_OTL_groups:
        dbDNA_df_OTL_groups.remove('')
    dbDNA_df_OTL = dbDNA_df.loc[dbDNA_df['species_group'].isin(dbDNA_df_OTL_groups)]
    dbDNA_df_OTL['OTL_taxon'] = dbDNA_df_OTL['species_name'].isin(OTL_species)
    dbDNA_OTL_xlsx = '/Volumes/Coruscant/dbDNA/MZB_Germany_v1_BarCodeBank/3_BarCodeBank/MZB_Germany_v1.BarCodeBank_OTL.xlsx'
    dbDNA_df_OTL.to_excel(dbDNA_OTL_xlsx, index=False)

########################################################################################################################

def main():

    # settings_xlsx = '/Volumes/Coruscant/dbDNA/settings/settings_mzb_mac.xlsx'

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
        data_source = tasks.loc[tasks['Task'] == 'source']['Run'].values.tolist()[0]
        run_download = tasks.loc[tasks['Task'] == 'download']['Run'].values.tolist()[0]
        run_extraction = tasks.loc[tasks['Task'] == 'extract']['Run'].values.tolist()[0]
        run_blacklist = tasks.loc[tasks['Task'] == 'blacklist']['Run'].values.tolist()[0]
        run_phylogeny = tasks.loc[tasks['Task'] == 'phylogeny']['Run'].values.tolist()[0]
        run_rating = tasks.loc[tasks['Task'] == 'rating']['Run'].values.tolist()[0]
        run_create_database = tasks.loc[tasks['Task'] == 'create database']['Run'].values.tolist()[0]
        run_create_report = tasks.loc[tasks['Task'] == 'create report']['Run'].values.tolist()[0]

        # Collect variables from settings file
        variables = pd.read_excel(settings_xlsx, sheet_name='Variables')
        project_name_or = variables.loc[variables['Variable'] == 'project name']['User input'].values.tolist()[0]
        version = variables.loc[variables['Variable'] == 'version']['User input'].values.tolist()[0]
        project_name = project_name_or + f'_v{version}'
        taxa_list = variables.loc[variables['Variable'] == 'taxa list']['User input'].values.tolist()[0]
        identifier_whitelist = variables.loc[variables['Variable'] == 'identifier whitelist']['User input'].values.tolist()[0]
        location_whitelist = variables.loc[variables['Variable'] == 'location whitelist']['User input'].values.tolist()[0]
        record_blacklist = variables.loc[variables['Variable'] == 'record blacklist']['User input'].values.tolist()[0]
        output_folder = variables.loc[variables['Variable'] == 'output folder']['User input'].values.tolist()[0]
        marker = variables.loc[variables['Variable'] == 'marker']['User input'].values.tolist()[0]
        mafft_executable = variables.loc[variables['Variable'] == 'mafft executable']['User input'].values.tolist()[0]
        similarity_threshold = variables.loc[variables['Variable'] == 'similarity threshold']['User input'].values.tolist()[0]
        vsearch_executable = variables.loc[variables['Variable'] == 'vsearch executable']['User input'].values.tolist()[0]
        makeblastdb_exe = variables.loc[variables['Variable'] == 'makeblastdb executable']['User input'].values.tolist()[0]
        midori2_fasta = variables.loc[variables['Variable'] == 'MIDORI2 fasta']['User input'].values.tolist()[0]
        lat_db = variables.loc[variables['Variable'] == 'lat']['User input'].values.tolist()[0]
        lon_db = variables.loc[variables['Variable'] == 'lon']['User input'].values.tolist()[0]
        d1 = variables.loc[variables['Variable'] == 'distance1']['User input'].values.tolist()[0]
        d2 = variables.loc[variables['Variable'] == 'distance2']['User input'].values.tolist()[0]
        d3 = variables.loc[variables['Variable'] == 'distance3']['User input'].values.tolist()[0]
        use_coordinates = variables.loc[variables['Variable'] == 'coordinates']['User input'].values.tolist()[0]
        use_country = variables.loc[variables['Variable'] == 'country']['User input'].values.tolist()[0]
        cpu_count = variables.loc[variables['Variable'] == 'cpu count']['User input'].values.tolist()[0]
        reference_xlsx = variables.loc[variables['Variable'] == 'reference_xlsx']['User input'].values.tolist()[0]
        reference_taxon = variables.loc[variables['Variable'] == 'reference_taxon']['User input'].values.tolist()[0]

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
        if type(cpu_count) != int or cpu_count == 0:
            cpu_count = multiprocessing.cpu_count() - 1
            print(f'{datetime.now().strftime("%H:%M:%S")} - Automatically detecing CPUs: Found {available_cores} (-1).')
        elif cpu_count > available_cores:
            cpu_count = multiprocessing.cpu_count() - 1
            print('{} - Not enough CPUs available. Defaulting to {} CPUs instead.'.format(datetime.now().strftime("%H:%M:%S"), cpu_count))
        ## print CPU usage
        print('{} - Multithreading mode: Using up to {} CPUs per individual command.'.format( datetime.now().strftime("%H:%M:%S"), cpu_count))

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

        if run_blacklist == 'Yes':
            blacklist_filter(output_directories, record_blacklist)

        if run_phylogeny == 'yes':
            phylogenetic_approach(output_directories, mafft_executable, vsearch_executable, similarity_threshold, cpu_count)

        if run_rating == 'yes':
            rating_system(output_directories, identifier_whitelist, location_whitelist, project_name, cpu_count)

        if run_create_database == 'yes':
            create_database(output_directories, project_name, makeblastdb_exe)

        if run_create_report == 'yes':
            create_report(output_directories, project_name, taxa_list)

        ## close the log file
        print('{} - Writing to log file...'.format(datetime.now().strftime("%H:%M:%S")))

        ## finish script
        print('\n{} - Done. Have a nice day!\n'.format(datetime.now().strftime("%H:%M:%S")))

if __name__ == "__main__":
    main()