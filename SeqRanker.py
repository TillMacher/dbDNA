import sys
import time, re
import numpy as np
from joblib import Parallel, delayed
import requests
from pathlib import Path
import pandas as pd
from tqdm import tqdm
import glob, os
import subprocess
from datetime import datetime
import multiprocessing
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import plotly.graph_objects as go
from Bio import Entrez
from requests_html import HTMLSession
from requests.packages.urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter
from ete3 import NCBITaxa
from ete3 import Tree
import xmltodict, pickle

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

## BOLD

def boldsystems_api(taxon_name, dir_out, marker, i, n_taxa):

    ## Retry if download failed
    max_retries = 15

    # remove subspecies
    taxon_name = ' '.join(taxon_name.split(' ')[:2])

    # Adjust query
    query_taxon = taxon_name.replace(" ", "%20")

    # Define output folder
    output = Path('{}/{}.tsv'.format(dir_out, taxon_name))

    # Check if file already exists
    if os.path.isfile(output):
        # Try to open file
        try:
            tmp = pd.read_table(output, low_memory=False)
            # Check if download was successful
            if tmp.values.tolist()[0][0] == 'MDB2 Error: connect failed':
                os.remove(output)
                run_api = True
                print('{} - Error in {} - Retry.'.format(datetime.now().strftime("%H:%M:%S"), taxon_name, i, n_taxa))
            elif tmp.columns.tolist()[0] != 'processid':
                os.remove(output)
                run_api = True
                print('{} - Error in {} - Retry.'.format(datetime.now().strftime("%H:%M:%S"), taxon_name, i, n_taxa))
            else:
                # Do not run if existing file looks good
                run_api = False
                print('{} - {} was already downloaded ({}/{}).'.format(datetime.now().strftime("%H:%M:%S"), taxon_name, i, n_taxa))
        except:
            # remove erroneous file and download again
            os.remove(output)
            run_api = True
    else:
        run_api = True

    if run_api == True:
        # Define the URL of the file
        url = "http://v4.boldsystems.org/index.php/API_Public/combined?taxon={}&format=tsv".format(query_taxon)

        for attempt in range(max_retries):
            try:
                # Send a GET request to the URL
                response = requests.get(url)

                # Check if the request was successful
                if response.status_code == 200:
                    if response.text != '':
                        # Open a file in write mode
                        with open(output, "w") as file:
                            # Write the content of the response to the file
                            file.write(response.text)
                    print('{} - Finished {} ({}/{}).'.format(datetime.now().strftime("%H:%M:%S"), taxon_name, i, n_taxa))
                    time.sleep(1)
                    break
                elif response.status_code == 429:
                    print('Too many requests, sleeping...')
                    time.sleep(2**(attempt + 1))  # exponential back-off
                    continue
            except Exception as e:
                if attempt < max_retries - 1:  # i.e. not the last attempt
                    time.sleep(2**(attempt + 1))  # exponential back-off
                    continue
                else:
                    print('{} - Failed to download {} after {} attempts.'.format(datetime.now().strftime("%H:%M:%S"), taxon_name, max_retries))
                    break

def download_data_from_bold(taxa_list, output_directories, marker):

    print('{} - Starting data download from bold.'.format(datetime.now().strftime("%H:%M:%S")))

    taxa_list_df = pd.read_excel(taxa_list).fillna('')
    dir_out = output_directories[1]

    taxa = [i[0] for i in taxa_list_df.values.tolist() if i[0] != '']
    n_taxa = len(taxa)

    # Download files via API
    Parallel(n_jobs = 5, backend='threading')(delayed(boldsystems_api)(taxon_name[1], dir_out, marker, taxon_name[0]+1, n_taxa) for taxon_name in enumerate(taxa))

    print('{} - Finished data download from bold.\n'.format(datetime.now().strftime("%H:%M:%S")))

def extract_bold_tsvs(output_directories, marker):

    print('{} - Starting data extraction from bold files.'.format(datetime.now().strftime("%H:%M:%S")))

    folder = output_directories[1]
    files = glob.glob('{}/*'.format(folder))
    all_record_list = []

    for file in files:
        df = pd.read_table(file, low_memory=False).fillna('')
        df = df.replace({' ': '', '': '', float('NaN'): ''})
        df = df.loc[df['species_name'] != ''] # species name is required
        df = df.loc[df['family_name'] != ''] # family name is required
        if len(df) != 0:
            all_record_list.extend(df.loc[df['nucleotides'] != ''].values.tolist())
        else:
            print('{} - Warning: No species data found for {} - species is removed.'.format(datetime.now().strftime("%H:%M:%S"), Path(file).stem))

    cols = pd.read_table(files[0]).fillna('').columns.tolist()
    df = pd.DataFrame(all_record_list, columns=cols)
    df_filtered = df.loc[df['markercode'] == marker] # the API still downloads other markers - so filter also at this step
    df_filtered = df_filtered.astype('string') # convert to strings, otherwise parquet crashes
    df_filtered = df_filtered.drop_duplicates() # drop potential duplicates
    output_file_1 = '{}/{}.parquet.snappy'.format(output_directories[3], 'file_1')
    df_filtered.to_parquet(output_file_1)

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

    taxids_xlsx = '/Users/tillmacher/Desktop/Projects/dbDNA/taxids.xlsx'
    taxids_df = pd.read_excel(taxids_xlsx).drop_duplicates()
    taxids_dict = {i[0]:i[1:] for i in taxids_df.values.tolist()}

    files = glob.glob(f'{output_directories[1]}/*.gb')

    columns = ['processid', 'sampleid', 'recordID', 'catalognum', 'fieldnum', 'institution_storing', 'collection_code', 'bin_uri',
     'phylum_taxID', 'phylum_name', 'class_taxID', 'class_name', 'order_taxID', 'order_name', 'family_taxID',
     'family_name', 'subfamily_taxID', 'subfamily_name', 'genus_taxID', 'genus_name', 'species_taxID', 'species_name',
     'subspecies_taxID', 'subspecies_name', 'identification_provided_by', 'identification_method',
     'identification_reference', 'tax_note', 'voucher_status', 'tissue_type', 'collection_event_id', 'collectors',
     'collectiondate_start', 'collectiondate_end', 'collectiontime', 'collection_note', 'site_code',
     'sampling_protocol', 'lifestage', 'sex', 'reproduction', 'habitat', 'associated_specimens', 'associated_taxa',
     'extrainfo', 'notes', 'lat', 'lon', 'coord_source', 'coord_accuracy', 'elev', 'depth', 'elev_accuracy',
     'depth_accuracy', 'country', 'province_state', 'region', 'sector', 'exactsite', 'image_ids', 'image_urls',
     'media_descriptors', 'captions', 'copyright_holders', 'copyright_years', 'copyright_licenses',
     'copyright_institutions', 'photographers', 'sequenceID', 'markercode', 'genbank_accession', 'nucleotides',
     'trace_ids', 'trace_names', 'trace_links', 'run_dates', 'sequencing_centers', 'directions', 'seq_primers',
     'marker_codes']

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
    df = pd.DataFrame(all_references_list, columns=columns)
    output_snappy = f'{output_directories[3]}/file_1.parquet.snappy'
    df.to_parquet(output_snappy, index=False)

    # update taxids file
    taxids_dict_2 = pd.DataFrame([[key] + values for key,values in taxids_dict.items()], columns=taxids_df.columns)
    taxids_dict_2.to_excel(taxids_xlsx, index=False)

    print('{} - Finished to collect data from .gb files.\n'.format(datetime.now().strftime("%H:%M:%S")))

## Phylogeny
# run alignment calculation
def run_mafft(df, family, mafft_executable, folder, outgroup_fasta):

    species_data = []
    sp = 0

    # family = 'Aeolosomatidae'
    # outgroup_fasta = '/Users/tillmacher/Downloads/fasta.fas'

    # Sub directories for output and temporary files
    family_dir = Path('{}/{}'.format(folder, family))
    if not os.path.exists(family_dir):
        os.makedirs(family_dir)

    ## required files to create
    fasta_file = '{}/{}.fasta'.format(family_dir, family)
    aln_file = '{}/{}.aln'.format(family_dir, family)
    aln_file_reduced = aln_file + '.reduced'
    aln_file_reduced_pickle = aln_file + '.reduced.pkl'
    species_file = '{}/{}.species'.format(family_dir, family)
    species_file_txt = species_file + ".txt"
    species_file_txt_snappy = Path(species_file_txt + ".parquet.snappy")

    f_out = open(f'{family_dir}/0_mafft.stdout.txt', 'w')
    f_err = open(f'{family_dir}/0_mafft.stderr.txt', 'w')

    # Check if file already exists
    if os.path.isfile(species_file_txt_snappy):
        message = '{} - {} was already analysed!'.format(datetime.now().strftime("%H:%M:%S"),family)

    else:
        ## Create fasta file
        bold_records = df.loc[df['family_name'] == family].values.tolist()

        # at least three sequences are required for the alignment!
        if len(bold_records) < 4:

            # create fasta file
            f = open(species_file_txt, 'w')
            for record in bold_records:
                f.write('{}__{}\n'.format(record[68], record[21]))
            f.close()
            message =  '{} - Skipped alignment for {} (less than 4 sequences).'.format(datetime.now().strftime("%H:%M:%S"),family)

        elif os.path.isfile(aln_file) and os.path.getsize(aln_file) > 0:
            message = '{} - Alignment was already calculated for {}.'.format(datetime.now().strftime("%H:%M:%S"),family)

        else:
            ############################################################################################################
            ## 1 - create fasta file

            # save outgroup names
            outgroup = []

            # create fasta file
            f = open(fasta_file, 'w')

            # enumerate to create unique identifier, which BOLD does not manage to do
            for i, record in enumerate(bold_records):

                # remove records with ambigous assignments
                species = record[21]
                special_characters = '!@#$%^&*()-+?_=,.<>/\'\"0123456789'
                # search for special characters in the species name
                if any(c in special_characters for c in species):
                    species_data.append([record[68], record[21], 0, 'Ambiguous assignment'])
                else:
                    ## format the fasta header
                    header = '>{}__{}__{}\n'.format(record[68], record[21].replace(" ", '_'), i)

                    ## remove all possible characters that infer with the newick format
                    characters_to_replace = [" ", ";", ":", ",", "(", ")", "'"]
                    for char in characters_to_replace:
                        header = header.replace(char, "")

                    ## write header
                    f.write(header)

                    ## write sequence
                    sequence = '{}\n'.format(record[71])
                    f.write(sequence)

            ## also add the sequence from the outgroup file
            for i, record in enumerate(SeqIO.parse(outgroup_fasta, "fasta")):
                name = re.sub('\W+', '', record.id)
                header = f">{name}_outgroup_{i}\n"

                ## remove all possible characters that infer with the newick format
                characters_to_replace = [" ", ";", ":", ",", "(", ")", "'"]
                for char in characters_to_replace:
                    header = header.replace(char, "")

                ## write header
                f.write(header)
                outgroup.append(f"{name}_outgroup_{i}")

                ## write sequence
                seq = f"{str(record.seq)}\n"
                f.write(seq)
            f.close()

            ############################################################################################################
            ## 2 - create alignemnt using mafft

            # create alignment
            command = "{} --auto --preservecase {} > {}".format(mafft_executable, fasta_file, aln_file)
            # Execute the command and capture stdout and stderr
            process = subprocess.Popen(command, shell=True, stdout=f_out, stderr=f_err)
            process.wait()

            message = '{} - Finished mafft alignment for {}.'.format(datetime.now().strftime("%H:%M:%S"), family)

        ############################################################################################################
        ## remove duplicates from alignment

        if os.path.isfile(aln_file) and os.path.getsize(aln_file) > 0 and not os.path.isfile(aln_file_reduced):
            ## reduce fasta to reduce iqtree2 runtimes
            reduced_fasta_dict = {}
            with open(aln_file) as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    id = record.id
                    seq = str(record.seq).upper()
                    if seq not in list(reduced_fasta_dict.keys()):
                        reduced_fasta_dict[seq] = [id]
                    else:
                        reduced_fasta_dict[seq] = reduced_fasta_dict[seq] + [id]

            ## write reduce fasta file
            ids_pickle = {}
            f = open(aln_file_reduced, 'w')
            i = 0
            for seq, ids in reduced_fasta_dict.items():
                ## replace header with temporary name
                tmp_header = f'>sp_{i}\n'
                sequence = f'{seq}\n'
                ## store true names in pickle file
                records = ';;'.join(ids)
                tmp_id = f'sp_{i}'
                ids_pickle[tmp_id] = records
                ## write to fasta
                f.write(tmp_header)
                f.write(sequence)
                ## increase counter
                i += 1
            f.close()

            with open(aln_file_reduced_pickle, 'wb') as handle:
                pickle.dump(ids_pickle, handle, protocol=pickle.HIGHEST_PROTOCOL)

            message = '{} - Finished (reduced) mafft alignment for {}.'.format(datetime.now().strftime("%H:%M:%S"),family)

    f_err.close()
    f_out.close()

    print(message)

## run phylogenetic tree calculation
def run_phylo(df, family, iqtree_executable, folder, cpu_count):

    # family = 'Chironomidae'
    # outgroup_fasta = '/Volumes/Coruscant/dbDNA/outgroup.fasta'

    # Sub directories for output and temporary files
    family_dir = Path('{}/{}'.format(folder, family))
    if not os.path.exists(family_dir):
        os.makedirs(family_dir)

    ## required files to create
    aln_file = '{}/{}.aln'.format(family_dir, family)
    aln_file_reduced = aln_file + '.reduced'
    tree_file = aln_file_reduced + '.treefile'
    species_file = '{}/{}.species'.format(family_dir, family)
    species_file_txt = species_file + ".txt"
    species_file_txt_snappy = Path(species_file_txt + ".parquet.snappy")

    ## collect bold records
    bold_records = df.loc[df['family_name'] == family].values.tolist()

    f_out = open(f'{family_dir}/1_phylo.stdout.txt', 'w')
    f_err = open(f'{family_dir}/1_phylo.stderr.txt', 'w')

    # Check if file already exists
    if os.path.isfile(species_file_txt_snappy):
        message = '{} - {} was already analysed!'.format(datetime.now().strftime("%H:%M:%S"), family)

    elif not os.path.isfile(aln_file_reduced):
        message = '{} - {} is missing an alignment file. Cannot calculate tree!'.format(datetime.now().strftime("%H:%M:%S"), family)

        # create fasta file
        f = open(species_file_txt, 'w')
        for record in bold_records:
            f.write('{}__{}\n'.format(record[68], record[21]))
        f.close()

    else:
        ## create tree using iqtree2 fast ##

        if os.path.isfile(tree_file):
            message = '{} - {} already has a tree file!'.format(datetime.now().strftime("%H:%M:%S"), family)

        else:
            ## calculate tree using iqtree2
            command = [iqtree_executable, '-s', aln_file_reduced, '-m', 'K2P', '-T', cpu_count, '--ninit', '1', '--fast']
            subprocess.call(command, stdout=f_out, stderr=f_err)

            message = '{} - Finished phylogenetic tree for {}.'.format(datetime.now().strftime("%H:%M:%S"),family)

    f_err.close()
    f_out.close()

    print(message)

## run species delnineation
def run_mptp(df, family, mptp_executable, folder, outgroup_fasta):

    # family = 'Branchinectidae'
    # outgroup_fasta = '/Users/tillmacher/Downloads/fasta.fas'

    species_data = []
    sp = 0

    # Sub directories for output and temporary files
    family_dir = Path('{}/{}'.format(folder, family))
    if not os.path.exists(family_dir):
        os.makedirs(family_dir)

    ## required files to create
    aln_file = '{}/{}.aln'.format(family_dir, family)
    aln_file_reduced = aln_file + '.reduced'
    aln_file_reduced_pickle = aln_file + '.reduced.pkl'
    tree_file = aln_file_reduced + '.treefile'
    species_file = '{}/{}.species'.format(family_dir, family)
    species_file_txt = species_file + ".txt"
    species_file_txt_snappy = Path(species_file_txt + ".parquet.snappy")

    ## collect bold records
    bold_records = df.loc[df['family_name'] == family].values.tolist()

    ## load ids and outgroup temporary names
    # Load pickel file
    with open(aln_file_reduced_pickle, 'rb') as handle:
        ids_dict = pickle.load(handle)
    outgroup = list(ids_dict.keys())[-1]

    ## open log files
    f_out = open(f'{family_dir}/2_mptp.stdout.txt', 'w')
    f_err = open(f'{family_dir}/2_mptp.stderr.txt', 'w')

    # Check if file already exists
    if os.path.isfile(species_file_txt_snappy):
        print('{} - {} was already analysed!'.format(datetime.now().strftime("%H:%M:%S"),family))

    elif not os.path.isfile(tree_file):
        message = '{} - {} is missing a tree file!'.format(datetime.now().strftime("%H:%M:%S"), family)

        # create fasta file
        f = open(species_file_txt, 'w')
        for record in bold_records:
            f.write('{}__{}\n'.format(record[68], record[21]))
        f.close()

    else:
        # identify species using mptp
        command = f"{mptp_executable} --ml --multi --tree_file {str(tree_file)} --outgroup {outgroup} --outgroup_crop --output {str(species_file)}"
        subprocess.run(command, shell=True, stdout=f_out, stderr=f_err)

        ## sleep to give time to write the file (happened several times that the file was not yet properly written)
        time.sleep(10)

        ## check if a species delineation file was created
        if not os.path.isfile(species_file_txt):
            # create fasta file
            f = open(species_file_txt, 'w')
            for record in bold_records:
                f.write('{}__{}\n'.format(record[68], record[21]))
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
                        tmp_id = line.strip()
                        ids = ids_dict[tmp_id].split(';;')

                        ## loop through all ids (can be more than one records, depending on the duplicates)
                        for record in ids:
                            ## collect original sequence information
                            species_id = record.split('__')[0]
                            species_name = record.split('__')[1]
                            species_data.append([species_id, species_name, current_species, LRT])

    ## check if a species delineation file was created
    if os.path.isfile(species_file_txt_snappy):
        ## write output file
        species_df = pd.DataFrame(species_data, columns=['species_id', 'species_name', 'current_species', 'LRT'])
        species_df = species_df.drop_duplicates()
        species_df = species_df.astype('string') # Convert to string, otherwise parquet crashes
        species_df.to_parquet(species_file_txt_snappy, index=False)
        message = '{} - Finished {}'.format(datetime.now().strftime("%H:%M:%S"), family)

    else:
        ## write a dummy file if mptp failed
        species_data = [[record[68], record[21], 0, 'Error'] for record in bold_records]
        species_df = pd.DataFrame(species_data, columns=['species_id', 'species_name', 'current_species', 'LRT'])
        species_df = species_df.astype('string') # Convert to string, otherwise parquet crashes
        species_df.to_parquet(species_file_txt_snappy, index=False)
        message = '{} - Finished {}, but failed species delineation'.format(datetime.now().strftime("%H:%M:%S"), family)

    f_out.close()
    f_err.close()

    print(message)

# main script for phylogenetic approach
def phylogenetic_approach(output_directories, mafft_executable, iqtree_executable, mptp_executable):

    print('{} - Starting phylogenetic approach.'.format(datetime.now().strftime("%H:%M:%S")))

    # Store fasta, alignments and tree
    folder = output_directories[2]

    # Load data
    output_file_1 = '{}/{}.parquet.snappy'.format(output_directories[3], 'file_1')
    df = pd.read_parquet(output_file_1).fillna('')

    # Extract families
    families = sorted(set(df['family_name'].values.tolist()))

    # for testing:
    # families = ['Chironomidae']
    # files = glob.glob('/Volumes/Coruscant/dbDNA/FEI_genera_BarCodeBank/2_phylogeny/*/*.species.txt.parquet.snappy')
    # [os.remove(file) for file in files if os.path.isfile(file)]

    # Use most CPUs
    cpu_count = multiprocessing.cpu_count() -1

    print('{} - Using {}/{} CPUs.\n'.format(datetime.now().strftime("%H:%M:%S"), cpu_count, cpu_count+1))

    ## NEW VERSION!
    # 1) create alignments using parallel and a single core per family
    Parallel(n_jobs = cpu_count, backend='threading')(delayed(run_mafft)(df, family, mafft_executable, folder, outgroup_fasta) for family in families)
    print('{} - Finished mafft alignments.\n'.format(datetime.now().strftime("%H:%M:%S")))

    # 2) calculate trees using all available cores per family - not in parallel
    [run_phylo(df, family, iqtree_executable, folder, str(cpu_count)) for family in families]
    print('{} - Finished calculation of phylogenetic trees.\n'.format(datetime.now().strftime("%H:%M:%S")))

    # 3) perform species delineation
    Parallel(n_jobs = cpu_count, backend='threading')(delayed(run_mptp)(df, family, mptp_executable, folder, outgroup_fasta) for family in families)
    print('{} - Finished species delineation.\n'.format(datetime.now().strftime("%H:%M:%S")))

    # Sleep one sec
    time.sleep(1)

    print('{} - Finished phylogenetic approach.\n'.format(datetime.now().strftime("%H:%M:%S")))

    print('{} - Starting to rate sequences (phylogeny).\n'.format(datetime.now().strftime("%H:%M:%S")))

    # Create a DataFrame
    lst = []
    for family in families:
        lst.extend(pd.read_parquet('{}/{}/{}.species.txt.parquet.snappy'.format(folder, family, family)).values.tolist())
    df_species = pd.DataFrame(lst, columns=['Sequence ID', 'Species Name', 'Species Group', 'LRT'])

    # Test Phylogeny
    res = []
    for species_group in set(df_species['Species Group'].values.tolist()):
        sub_df = df_species.loc[df_species['Species Group'] == species_group]
        taxon = list(set(sub_df['Species Name'].values.tolist()))
        LRT = list(set(sub_df['LRT'].values.tolist()))[0]
        sub_df = sub_df.drop('LRT', axis=1)
        n_records = len(sub_df)

        ## monophyletic and more than 2 sequences == monophyletic (1)
        if LRT == 'Passed' and n_records >= 2:
            phylogeny = 'monophyletic'

        ## monophyletic but single sequence == monophyletic (2)
        elif LRT == 'Passed' and n_records < 2:
            phylogeny = 'monophyletic (singleton)'

        ## LRT failed == paraphyletic (1)
        elif LRT == 'Failed':
            phylogeny = 'paraphyletic (LRT failed)'

        ## LRT insufficient data == paraphyletic (2)
        else:
            phylogeny = 'paraphyletic (insufficient data)'

        ## add to results
        for record in sub_df.values.tolist():
            res.append(record + [0, phylogeny, 0])

    df_out = pd.DataFrame(res, columns=['Sequence ID', 'Species Name', 'Species Group', 'Clade', 'State', 'Rating'])
    df_out = df_out.astype('string')  # Convert to string, otherwise parquet crashes
    output_file_2 = '{}/{}.parquet.snappy'.format(output_directories[3], 'file_2')
    df_out.to_parquet(output_file_2, index=False)

    print('{} - Finished to rate sequences (phylogeny).\n'.format(datetime.now().strftime("%H:%M:%S")))

## Rating system

def rating_system(output_directories, identifier_whitelist, location_whitelist, project_name):

    print('{} - Starting to rate sequences.'.format(datetime.now().strftime("%H:%M:%S")))

    all_ratings_list = []

    # Identifier white list
    identifier_white_list = pd.read_excel(identifier_whitelist).fillna('')['identifier_white_list'].values.tolist()

    # Location white list
    country_white_list = pd.read_excel(location_whitelist).fillna('')
    main_country = [i for i in country_white_list['main country'].values.tolist() if i != '']
    neighbour_countries = [i for i in country_white_list['neighbour_countries'].values.tolist() if i != '']
    continent = [i for i in country_white_list['continent'].values.tolist() if i != '']

    # Collect phylogenetic approach values
    file_1 = Path('{}/file_1.parquet.snappy'.format(output_directories[3]))
    df1 = pd.read_parquet(file_1).fillna('')
    file_2 = Path('{}/file_2.parquet.snappy'.format(output_directories[3]))
    df2 = pd.read_parquet(file_2).fillna('')

    t = 0
    for record in tqdm(df2.values.tolist()):
        t += 1
        if t == 15000:
            break

        # Collect information for the specific record

        # Phylogeny file
        sequenceID = record[0]
        species_group = record[2]
        clade = record[3]
        phylogeny = record[4]
        rating = int(record[5])

        # BOLD record file
        bold_record = df1.loc[df1['sequenceID'] == sequenceID].values.tolist()[0]
        processid = bold_record[0]
        institution_storing = bold_record[5]
        bin_uri = bold_record[7]
        phylum_name = bold_record[9]
        class_name = bold_record[11]
        order_name = bold_record[13]
        family_name = bold_record[15]
        genus_name = bold_record[19]
        species_name = bold_record[21]
        identification_by = bold_record[24]
        country = bold_record[54]
        province = bold_record[55]
        region = bold_record[56]
        exactsite = bold_record[58]
        lifestage = bold_record[38]
        sex = bold_record[39]
        image_urls = bold_record[60]
        markercode = bold_record[69]
        nucleotides = bold_record[71]
        sequence_quality = 'Sufficient'

        ## phylogeny
        if phylogeny == 'monophyletic':
            rating += 15
        elif phylogeny == 'monophyletic (singleton)':
            rating += 5

        ## trimm leading and trailing gaps from barcode
        nucleotides_trimmed = nucleotides.strip('-').strip('N')

        ## good sequence quality
        allowed_chars = set('ACGT')
        if set(nucleotides_trimmed).issubset(allowed_chars):
            rating += 5
            sequence_quality = 'Good'


        ## bad sequence quality
        not_allowed_chars = [i for i in set(nucleotides_trimmed) if i not in allowed_chars]
        if len(not_allowed_chars) != 0:
            n_not_allowed_chars = sum([nucleotides_trimmed.count(i) for i in not_allowed_chars])
            rel = n_not_allowed_chars / len(nucleotides_trimmed) * 100
            if rel >= 2:
                rating -= 10
                sequence_quality = 'Bad'

        ## sequence length (>= 500 bp are accepted as barcode)
        if len(nucleotides_trimmed) >= 500:
            rating += 5
            barcode_length = len(nucleotides_trimmed)

        ## Identification white list
        if identification_by in identifier_white_list:
            rating += 10

        ## Sampling location with white list
        if country in main_country:
            rating += 5
        elif country in neighbour_countries:
            rating += 4
        elif country in continent:
            rating += 3

        ## Check if photo is available
        if image_urls != '':
            rating += 5

        ## Available metadata
        if province != '':
            rating += 1
        if region != '':
            rating += 1
        if exactsite != '':
            rating += 1
        if lifestage != '':
            rating += 1
        if sex != '':
            rating += 1

        all_ratings_list.append([rating, sequenceID, processid, bin_uri, phylum_name, class_name, order_name,
                                 family_name, genus_name, species_name, clade, phylogeny, species_group, identification_by,
                                 institution_storing, country, province, region, exactsite, lifestage, sex, image_urls, markercode,
                                 sequence_quality, barcode_length, nucleotides])

    # create unfiltered dataframe
    ratings_df = pd.DataFrame(all_ratings_list, columns=["rating", "sequenceID", "processid", "bin_uri", "phylum_name", "class_name", "order_name",
                                 "family_name", "genus_name", "species_name", "clade", "phylogeny", "species_group", "identification_by",
                                 "institution_storing", "country", "province", "region", "exactsite", "lifestage", "sex", "image_urls", "markercode",
                                 "sequence_quality", "barcode_length", "nucleotides"])
    ratings_df = ratings_df.sort_values('rating', ascending=False)
    ratings_df = ratings_df.astype('string')  # Convert to string, otherwise parquet crashes
    ratings_df['rating'] = ratings_df['rating'].astype('int')
    ratings_df['clade'] = ratings_df['clade'].astype('int')

    # filter dataframe according to the rating threshold
    ratings_df_filtered = ratings_df.loc[ratings_df['rating'] >= 5]

    # write to parquet
    output_file_1 = Path('{}/{}.BarCodeBank.parquet.snappy'.format(output_directories[3], project_name))
    ratings_df.to_parquet(output_file_1, index=False)
    output_file_2 = Path('{}/{}.BarCodeBank_filtered.parquet.snappy'.format(output_directories[3], project_name))
    ratings_df_filtered.to_parquet(output_file_2, index=False)

    # try to write to excel
    try:
        output_file_1 = Path('{}/{}.BarCodeBank.xlsx'.format(output_directories[3], project_name))
        ratings_df.to_excel(output_file_1, index=False)
        output_file_2 = Path('{}/{}.BarCodeBank_filtered.xlsx'.format(output_directories[3], project_name))
        ratings_df_filtered.to_excel(output_file_2, index=False)
    except:
        print('{} - Unable to write to xlsx.'.format(datetime.now().strftime("%H:%M:%S")))

    print('{} - Finished to rate sequences.\n'.format(datetime.now().strftime("%H:%M:%S")))

def create_report():
    ratings_snappy = Path('{}/{}.BarCodeBank.parquet.snappy'.format(output_directories[3], project_name))
    ratings_df = pd.read_parquet(ratings_snappy)

    ratings_filtered_snappy = Path('{}/{}.BarCodeBank_filtered.parquet.snappy'.format(output_directories[3], project_name))
    ratings_df_filtered = pd.read_parquet(ratings_filtered_snappy)

    ##### PLOT 1 #####
    ## Rating distribution
    all_ratings = ratings_df['rating'].values.tolist()
    max_rating = 35 # ADJUST THIS IF THE POINTS CHANGE!!
    y = [all_ratings.count(i) for i in range(1,max_rating+1)]

    fig = go.Figure()
    # BRONZE
    fig.add_vrect(x0=4.5, x1=14.5, line_width=0, fillcolor="Peru", opacity=0.3, layer='below')
    # SILVER
    fig.add_vrect(x0=14.5, x1=24.5, line_width=0, fillcolor="Silver", opacity=0.3, layer='below')
    # GOLD
    fig.add_vrect(x0=24.5, x1=35.5, line_width=0, fillcolor="Gold", opacity=0.3, layer='below')
    fig.add_trace(go.Bar(x=[i for i in range(0,max_rating+1)], y=y, marker_color='navy', marker_line_color='navy'))
    fig.update_layout(template='simple_white', width=600, height=500)
    fig.update_yaxes(title='Reference sequences', dtick='linear')
    fig.update_xaxes(title='Rating', dtick='linear')
    fig.write_image('/Users/tillmacher/Desktop/Projects/dbDNA/Präsentationen/report_a.pdf')

    ##### PLOT 2 #####
    ## Number of sequences per species
    all_species = [i for i in ratings_df['species_name'].values.tolist() if i != 'No Match']
    res = {i:all_species.count(i) for i in set(all_species)}
    res = dict(sorted(res.items(), key=lambda item: item[1]))

    x = list(res.values())[-30:]
    y = list(res.keys())[-30:]

    fig = go.Figure()
    fig.add_trace(go.Bar(x=x, y=y, orientation='h'))
    fig.update_layout(template='simple_white', width=600, height=500)
    fig.update_xaxes(title='Reference sequences')
    fig.update_yaxes(dtick='linear')
    fig.write_image('/Users/tillmacher/Desktop/Projects/dbDNA/Präsentationen/report_b.pdf')

    ##### PLOT 3 #####
    ## Average rating per species
    all_species = ratings_df.sort_values('family_name', ascending=False)['species_name'].drop_duplicates().values.tolist()
    fig = go.Figure()
    res = {}
    for species in all_species:
        ratings = ratings_df.loc[ratings_df['species_name'] == species]['rating'].values.tolist()
        avg = np.average(ratings)
        res[species] = [avg, ratings]
    # sort the species to only show the x top species
    res = dict(sorted(res.items(), key=lambda item: item[1][0]))
    top_species = list(res.keys())[-50:]

    for species in top_species:
        ratings = res[species][1]
        fig.add_trace(go.Box(x=ratings, name=species, orientation='h'))

    fig.update_layout(template='simple_white', width=600, height=500)
    fig.update_xaxes(title='Reference sequences', range=(0,35))
    fig.write_image('/Users/tillmacher/Desktop/Projects/dbDNA/Präsentationen/report_c.pdf')

## BLASTn

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
    if os.path.isfile(f'{db_name}.ndb'):
        print('{}: Finished building database.'.format(datetime.now().strftime('%H:%M:%S')))
    else:
        print('{}: An error occurred when building the database.'.format(datetime.now().strftime('%H:%M:%S')))

########################################################################################################################

if len(sys.argv) > 1:
    settings_xlsx = Path(sys.argv[1])
else:
    settings_xlsx = Path('./settings.xlsx')

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

########################################################################################################################

new_project = 'yes'
if new_project == 'yes':
    output_directories = create_new_project(project_name, output_folder)

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
        extract_bold_tsvs(output_directories, marker)

    # GenBank workflow
    elif data_source == 'NCBI':
        extract_genbank_files(output_directories)

if run_phylogeny == 'yes':
    phylogenetic_approach(output_directories, mafft_executable, iqtree_executable, mptp_executable)

if run_rating == 'yes':
    rating_system(output_directories, identifier_whitelist, location_whitelist, project_name)

if run_create_database == 'yes':
    create_database(output_directories, project_name, makeblastdb_exe)



'./decenttree -fasta /Volumes/Coruscant/dbDNA/test/Limoniidae.aln -out /Volumes/Coruscant/dbDNA/test/Limoniidae.newick -t BIONJ'

'/Volumes/Coruscant/dbDNA/programs/mptp --ml --multi --tree_file /Volumes/Coruscant/dbDNA/test/Ameletidae.newick --output /Volumes/Coruscant/dbDNA/test/Ameletidae.species'
