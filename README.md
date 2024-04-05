# AMRDB - ResFinder Results stored in Database

This repository contains the python package for AMRDB with two main functionalities:
- Updating the internal sequence database based on ResFinder_db git repository
- Adding new results from several tools (ResFinder, Bakta, ISEScan2, AMRFinderPlus, SpeciesFinder,
 PlasmidFinder)
  
The general concept is, that detected resistance gene sequences are stored in a database format. Among those
 sequences, that were originally present in ResFinder-DB and AMRFinderPlus database (nomenclature-databases),
 the novel variants of each gene is stored. E.g. when resistance genes are detected but have less than 100%
 identity to the closest reference, this exact sequence receives a local variant-name, which will be
 complemented by the real name as soon as the sequence is added to the nomenclature-databases.  
  
Additional results from annotation tools like ISEScan2, Bakta and PlasmidFinder are supposed to provide
 insights into the transmissability and mobility of the contig/element the resistance conveying genes are
 located on.

A nextflow workflow integrating this database package can be found [here](https://github.com/ages-bioinformatics/wf-antibiotic-resistance)  


## Requirements

- Python3 (developed in 3.10.9)
- sqlalchemy>=2.0 (developed in 2.0.18)
- pandas (developed in 1.5.3)
- mysql-connector-python (developed in 8.0.33)

## Usage

For updating or initializing a database, run
`update_resfinder_database.py --resfinder_db /path/to/resfinder_db`

In default behaviour, this will create a sqlite-database `agres.db` in current workingdirectory. Database-path can be changed with `-d /path/to/database.db`

To connect a mysql or mariadb-database, specify `-d database_name`, `-H hostname`, `-u db_username` and `-p db_user_password`.

Results from tools can be added using `amrdb_add_result.py` with same database-parameters as above.  
Additionally, an input-path must be specified with `-i`  
Depending on the `--method` argument these are:


| `--method` |            `-i`              |
|------------|------------------------------|
| resfinder  | path/to/resfinder_output_dir |
| bakta      | path/to/bakta_result.tsv     |
| isescan    | path/to/isescan_result.tsv   |
| mobtyper   | path/to/mobtyper_result.tsv  |
| plasmidfinder  | path/to/plasmidfinder_results.tsv  |
| speciesfinder | path/to/data.txt          |
| amrfinder  | path/to/amrfinder_directory  |
| mlst       | path/to/mlst.csv             |

For ResFinder and AMRFinder, there is an option to import the fasta-file (`--assembly`). This will initialize and create contig elements in the database.  
Also, ResFinder output could be created from fastq-input. In this case use the `--mode fastq` parameter. This way it is not possible to extract coordinates of the found genes and `--assembly` input will not work.  
Note that method amrfinder takes a directory as input and expects two files named `amrfinder_results.txt`
and `amrfinder_nucleotides.fasta` in this directory.

## Usage of 3rd party tools

Bakta should be run using `--keep-contig-headers` flag.  
mob_suite input should be generated using `mob_typer`.  
AMRFinderPlus should be run with `--nucleotide_output amrfinder_output/amrfinder_nucleotides.fasta` and 
`-o amrfinder_output/amrfinder_results.txt` (and `--organism` if applicable), while `amrfinder_output` is the
directory used as input path for the import script.  
mlst should be run using `--csv` flag.  
Other tools can be run with default/user prefered parameters.  
