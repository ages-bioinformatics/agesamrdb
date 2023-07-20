#!/bin/bash

if [ ! -d resfinder_db ]; then
	git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git
fi

rm -f ./agres.db

python update_resfinder_database.py -H "$hostname" -d "$dbname" -u "$dbuser" -p "$dbpw" --resfinder_db "resfinder_db"

python add_result.py --method resfinder -i "testdata/resfinder_assembly_out" --external_id 1 --assembly testdata/NRLAR-22-ESBL-22-0863-1-FH20.fasta
python add_result.py --method resfinder -i "testdata/resfinder_reads_out" --external_id 1 --mode fastq
python add_result.py --method isescan -i "testdata/NRLAR-22-ESBL-22-0863-1-FH20.fasta.tsv" --external_id 1
python add_result.py --method bakta -i "testdata/bakta/NRLAR-22-ESBL-22-0863-1-FH20.tsv" --external_id 1
