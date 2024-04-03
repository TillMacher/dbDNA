# dbDNA - A phylogeny- and expert identifier-driven grading system for reliable taxonomic annotation of (meta)barcoding data

<img src="https://github.com/TillMacher/apscale_gui/blob/master/_data/apscale_start.png" width="40%" height="40%">

## Introduction

Text

## Installation
### SeqRanker
Individual dbDNA databases can be created using the SeqRanker pipeline, which can be installed on all common operating systems (Windows, Linux, MacOS). SeqRanker requires Python 3.7 or higher and can be easily installed via pip in any command line:

`pip3 install seqranker`

To update SeqRanker run:

`pip3 install --upgrade seqranker`

Besides the main script, several other programs are required for the database creation. Please follow the installation instructions for your operating system for each software.

### mafft
Mafft is software to calculate multiple sequence alignments and is required the phylogenetic approach. More information about the installation of mafft can be found [here](https://mafft.cbrc.jp/alignment/software/).

### IQ-TREE
IQ-TREE is a phylogenomic software that calculate maximum likelihood trees. IQ-TREE is required to for the phylogenetic approach. More information about the installation of IQ-TREE can be found [here](https://github.com/iqtree/iqtree2).

### mPTP
mPTP is a software that is applied for species delimitation using the multi-rate Poisson Tree Processes. More information about the installation of mPTP can be found [here](https://github.com/Pas-Kapli/mptp)

### BLAST+
BLAST+ is a software to create BLAST databases and perform BLAST searches on custom (local) databases. More information about the installation of BLAST+ can be found [here](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#downloadblastdata).

## Settings file

The SeqRanker pipeline collects the required information from a .xlsx file. All specifications must be entered into this file.

Sheet 1 contains the Run parameters. **Here, the "Run" column is to be modified**

| **Task**          | **Run** | **Comment**                       |
|:------------------|:-------:|:----------------------------------|
| source            | **BOLD**    | define source                     |
| download          | **yes**      | download BOLD/NCBI data          |
| extract           | **yes**      | extract BOLD/NCBI data           |
| phylogeny         | **yes**     | calculate phylogenetic trees      |
| rating            | **yes**     | create table and rate records     |
| create database   | **yes**     | create blast database             |

Sheet 2 contains the database information and source files. **Here, the "User input" column is to be modified**

| **Variable**          | **User input**          | **Comment**                               | **Options** |
|:----------------------|:------------------------|:------------------------------------------|:------------|
| project name          | **Invertebrate_example_database** | Name of the database                  | string      |
| taxa list             | **/PATH/invertebrates.xlsx** | Excel file containing taxa to download    | PATH        |
| identifier whitelist  | **/PATH/identifier_white_list.xlsx** | Enter path to identifier whitelist | PATH        |
| location whitelist    | **/PATH/country_white_list.xlsx** | Enter path to location whitelist   | PATH        |
| output folder         | **/PATH/example**           | Enter path to output directory             | PATH        |
| marker                | **COI-5P**                   | Marker to download                         | string      |
| rating minimum        | **5**                       | Keep only sequences that are >= X          | yes / no    |
| download overwrite    | **yes**                     | Overwrite existing files?                  | yes / no    |
| alignment overwrite   | **yes**                     | Overwrite existing files?                  | yes / no    |
| tree overwrite        | **yes**                     | Overwrite existing files?                  | yes / no    |
| mafft executable      | **/PATH/mafft**             | Either "mafft" or "PATH/TO/mafft"           | PATH        |
| iqtree executable     | **/PATH/iqtree2**           | Either "iqtree" or "PATH/TO/iqtree"         | PATH        |
| mptp executable       | **/PATH/mptp**              | Either "mptp" or "PATH/TO/mptp"             | PATH        |
| makeblastdb executable | **/PATH/makeblastdb**      | Either "makeblastdb" or "PATH/TO/makeblastdb" | PATH        |
| MIDORI2 fasta         |                          | Enter path to MDORI2 file                  | PATH        |
| outgroup_fasta        | **/PATH/outgroup.fasta**    | Enter path to outgroup sequence             | PATH        |


## Documentation

Text

## Available databases

Text

## Citation

Text
