#!/bin/bash

if [ ! -d resfinder_db ]; then
	git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git
fi

rm -f ./agres.db

python update_resfinder_database.py

python amrdb_add_results.py --method resfinder -i "testdata/resfinder_assembly_out" --external_id 1 --assembly testdata/NRLAR-22-ESBL-22-0863-1-FH20.fasta --sample_name test
python amrdb_add_results.py --method resfinder -i "testdata/resfinder_reads_out" --external_id 1 --mode fastq
python amrdb_add_results.py --method isescan -i "testdata/isescan/NRLAR-22-ESBL-22-0863-1-FH20.fasta.tsv" --external_id 1
python amrdb_add_results.py --method bakta -i "testdata/bakta/NRLAR-22-ESBL-22-0863-1-FH20.tsv" --external_id 1
python amrdb_add_results.py --method mobtyper -i "testdata/mobtyper_results.txt" --external_id 1
