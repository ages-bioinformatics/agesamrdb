# AMRDB - ResFinder Results stored in Database

This repository contains the python package for AMRDB with two main functionalities:
- Updating the internal sequence database based on ResFinder_db git repository
- Adding new results from several tools (ResFinder, Bakta, ISEScan2)

## Requirements

- Python3 (developed in 3.10.9)
- sqlalchemy>=2.0 (developed in 2.0.18)
- pandas (developed in 1.5.3) - version not critical
- mysql-connector-python (developed in 8.0.33)

## Usage

For updating or initializing a database, run
`update_resfinder_database.py --resfinder_db /path/to/resfinder_db`

In default behaviour, this will create a sqlite-database `agres.db` in current workingdirectory. Database-path can be changed with `-d /path/to/database.db`

To connect a mysql or mariadb-database, specify `-d database_name`, `-H hostname`, `-u db_username` and `-p db_user_password`.

Results from tools can be added using `add_result.py` with same database-parameters as above.  
Additionally, an input-path must be specified with `-i`  
Depending on the `--method` argument these are:


| `--method` |            `-i`              |
|------------|------------------------------|
| resfinder  | path/to/resfinder_output_dir |
| bakta      | path/to/bakta_result.tsv     |
| isescan    | path/to/isescan_result.tsv   |
| mobtyper   | path/to/mobtyper_result.tsv  |
| plasmidfinder  | path/to/plasmidfinder_results.tsv  |
| phispy     | path/to/prophage_coordinates.tsv  |

For ResFinder, there is an option to import the fasta-file (`--assembly`). This will initialize and create contig elements in the database. It is required to make best use of the visualisation webUI with dna_features_viewer.  
Also, ResFinder output could be created from fastq-input. In this case use the `--mode fastq` parameter. This way it is not possible to extract coordinates of the found genes and `--assembly` input will not work.

## Usage of 3rd party tools

Bakta should be run using `--keep-contig-headers` flag.  
mob_suite input should be generated using `mob_typer`.
Other tools can be run with default/user prefered parameters.  
A nextflow workflow integrating this database package can be found [here](to_be_added.html)  
