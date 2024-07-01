# dbDNA - A phylogeny- and expert identifier-driven grading system for reliable taxonomic annotation of (meta)barcoding data

## Introduction

Text

## Installation
#### SeqRanker pipeline
Individual dbDNA databases can be created using the SeqRanker pipeline, which can be installed on all common operating systems (Windows, Linux, MacOS). SeqRanker requires Python 3.7 or higher and can be easily installed via pip in any command line:

`pip3 install seqranker`

To update SeqRanker run:

`pip3 install --upgrade seqranker`

Alternatively, standalone versions of the SeqRanker pipeline for Windows11 and MacOS (tested on Ventura 13.5) are available under the [latest release](https://github.com/TillMacher/dbDNA/releases).

#### Further Dependencies

Besides the main script, several other programs are required for the database creation. Please follow the installation instructions for your operating system for each software.

#### mafft
Mafft is software to calculate multiple sequence alignments and is required the phylogenetic approach. More information about the installation of mafft can be found [here](https://mafft.cbrc.jp/alignment/software/).

#### IQ-TREE
IQ-TREE is a phylogenomic software that calculate maximum likelihood trees. IQ-TREE is required to for the phylogenetic approach. More information about the installation of IQ-TREE can be found [here](https://github.com/iqtree/iqtree2).

#### mPTP
mPTP is a software that is applied for species delimitation using the multi-rate Poisson Tree Processes. More information about the installation of mPTP can be found [here](https://github.com/Pas-Kapli/mptp)

#### BLAST+
BLAST+ is a software to create BLAST databases and perform BLAST searches on custom (local) databases. More information about the installation of BLAST+ can be found [here](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#downloadblastdata).

#### APSCALE blast
APSCALE is a software to process (e)DNA metabarcoding datasets. The blastn module is used to perform BLAST searches on custom (local) databases. More information about the installation of APSCALE blast can be found [here](https://github.com/TillMacher/apscale_gui).

## Settings file

The SeqRanker pipeline collects the required information from an [Excel file](https://github.com/TillMacher/dbDNA/blob/main/european_freshwater_invertebrates/settings.xlsx). All specifications must be entered into this file.

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

## Run SeqRanker
First, prepare the settings file according to your needs. Then, the SeqRanker pipeline can easily be initiated via the following command(s):

#### pypi version
* Open a new terminal
* Execute: `seqranker ./PATH/TO/FOLDER/settings.xlsx`

#### standalone version
* Doubleclick on the `seqranker_v0.1-macosx-ventura` or `seqranker_v0.1-W11` executable.
* Provide the settings.xlsx file.

## Example data
Example data that was used for the creation a database for European freshwater invertebrates can be found [here](https://github.com/TillMacher/dbDNA/tree/main/european_freshwater_invertebrates):
* [Taxa list](https://github.com/TillMacher/dbDNA/blob/main/european_freshwater_invertebrates/Freshwaterecology_info_all_invertebrates.xlsx)
* [Country white list](https://github.com/TillMacher/dbDNA/blob/main/european_freshwater_invertebrates/country_white_list.xlsx)
* [Identifier white list](https://github.com/TillMacher/dbDNA/blob/main/european_freshwater_invertebrates/identifier_white_list.xlsx)
* [Settings file](https://github.com/TillMacher/dbDNA/blob/main/european_freshwater_invertebrates/settings.xlsx)

## SeqRanker pipeline: a short overview

#### Overview slides
* A more detailed overview into the pipeline can be found in [this](https://github.com/TillMacher/dbDNA/blob/main/source/dbDNA_overview.pdf) presentation.

#### Step 1: Data acquisition
* Records for all taxa provided in taxa list are downloaded (the taxon can be any taxonomic level). For example, of a genus is provided, all species records for this genus will be fetched.
* Sequence records can be obtained from **BOLDsystems** and **MIDORI2** (GenBank).
* For each record, all available metadata is downloaded (from BOLDsystems or GenBank, depending on the source).
* All records and their respective metadata are stored in a raw sequence table.

#### Step 2: Species delineation
* The sequences of all records of each family in the dataset are combined in a separate .fasta file.
* A multiple sequence alignment for each family is calculated, using _mafft_.
* A maximum likelihood tree for each family is calculated, using _IQ-Tree_ (fast option).
* Species are delimited for each family, using _mPTP_.
* The species delimitation results are used evaluate if a species record is mono- or paraphyletic.

#### Step 3: Rating system
* Each individual record is scored, based on the following criteria.
* If a criterion is not met, no points are gained.

| **Category**          | **Points gained** | **Explanation**                               |
|:----------------------|:------------------:|:----------------------------------------------|
| monophyletic OR       | 15                | Delimited species group only contains one species |
| monophyletic (singleton) | 5               | Delimited species group only contains one species, but only a single sequence |
| good sequence quality | 3                 | Only the four bases "AGCT" are present         |
| bad sequence quality  | -10               | More than 2% of the sequence are not "AGCT"    |
| longer than 500 bp    | 2                 | The recommended minimum barcode length is >= 500 bp |
| identifier on whitelist | 15           | The specimen was identified by an identifier on the white list |
| main country OR       | 9                 | The specimen was collected in the main country |
| neighbour country OR  | 6                 | The specimen was collected in a neighbouring country |
| continent             | 3                 | The specimen was collected on the same continent |
| distance <= d1 OR     | 9                 | The specimen was collected in the main country |
| distance <= d2 OR     | 6                 | The specimen was collected in a neighbouring country |
| distance <= d3        | 3                 | The specimen was collected on the same continent |
| image                 | 1                 | An image is available                         |
| province              | 1                 | The metadata is available                            |
| region                | 1                 | The metadata is available                            |
| exactsite             | 1                 | The metadata is available                            |
| lifestage             | 1                 | The metadata is available                            |
| sex                   | 1                 | The metadata is available                            |

* Each record can gain between 50 (excellent) and -10 (highly unreliable) points.
* All records are categorized according to their points.

| **Border** | **Gold** | **Silver** | **Bronze** | **Unreliable** |
| --- | --- | --- | --- | --- |
| Upper | 50 | 39 | 24 | 9 |
| Lower | 40 | 25 | 10 | -10 |

#### Step 4: Database creation
* The function makeblastdb is used to create a BLAST+ compatible database.

#### Step 5: Local BLASTn
* The APSCALE BLASTn tool can be used for the taxonomic assignment of DNA metabarcoding datasets against the newly created database.
* APSCALE will automatically filter the hits and include the ratings of the record in the filtering process.
* The filtering algorithm works as follows, for each OTU individually:
1. Obtain the Top20 BLASTn hits for the OTU.
2. Filter by similarity: all hits with the highest similarity are kept.
3. Trim hits according to similarity: Species >=98%, Genus >=95%, Family >=90%, Order >= 85%.
4. Filter remaining hits by rating: A) keep all Gold hits OR B) keep all Silver hits OR C) keep all Bronze hits OR D) keep all unreliable hits.
5. Trim taxonomy of remaining hits to their most recent common ancestor (MRCA filtering): Phylum, Class, Order, Family, Genus, Species.
* All ambiguous taxonomic assignments and metadata are kept in the final table as "traits" for each OTU.

## Available databases
#### European freshwater invertebrates (COI)
* All species of all genera classified as European freshwater invertebrates (according to [freshwaterecology.info](https://www.freshwaterecology.info/index.php)).
* A filtered and unfilitered version is available [here](https://www.freshwaterecology.info/index.php).
#### European freshwater fish and lamprey (12S)
* All species of all genera classified as European freshwater fish and lamprey (according to [freshwaterecology.info](https://www.freshwaterecology.info/index.php)).
* A filtered and unfilitered version is available [here](https://www.freshwaterecology.info/index.php).

## Benchmark
* Runtimes for the SeqRanker database creation are optimized for parallelization.
* Increasing the number of available cores will signficantly reduce runtimes.
* However, even large databases can be curated on average hardware.

#### Example
* All genera of all European freshwater macroinvertebrates, available on freshwater-ecology.info.
* In total 500,521k records were downloaded from BOLDsystems.
* Executed on a MacBook M1 Pro 2021 (16GB RAM, 8 cores).

| **Runtime (min)** | **Step** |
| --- | --- |
| 124 | Sequence download |
| 2 | Record extraction |
| 20 | Alignments |
| 120 | ML tree |
| 10 | Species delimitation |
| 8 | Barcode ranking |
| 6 | Database creation |


## Citation

#### SeqRanker
_Coming soon..._

#### mafft
Katoh, K., Misawa, K., Kuma, K., & Miyata, T. (2002). MAFFT: A novel method for rapid multiple sequence alignment based on fast Fourier transform. Nucleic Acids Research, 30(14), 3059–3066. https://doi.org/10.1093/nar/gkf436

#### IQ-Tree
Nguyen, L.-T., Schmidt, H. A., von Haeseler, A., & Minh, B. Q. (2015). IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies. Molecular Biology and Evolution, 32(1), 268–274. https://doi.org/10.1093/molbev/msu300

#### mPTP
Kapli, P., Lutteropp, S., Zhang, J., Kobert, K., Pavlidis, P., Stamatakis, A., & Flouri, T. (2017). Multi-rate Poisson tree processes for single-locus species delimitation under maximum likelihood and Markov chain Monte Carlo. Bioinformatics, 33(11), 1630–1638. https://doi.org/10.1093/bioinformatics/btx025
