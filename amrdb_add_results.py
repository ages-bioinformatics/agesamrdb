#!/usr/bin/env python

import os
import argparse

from amrdb.models import Base, Sample 
from amrdb.interfaces import read_result, insert_into_db
from amrdb.util import get_or_create

from sqlalchemy import create_engine, inspect 
from sqlalchemy.orm import Session

# CLI definitions
parser = argparse.ArgumentParser(description="agres - helperscript to add " \
        + "single result to database - methods allowed are: " \
        + "resfinder, bakta, isescan, mobtyper")
# db definitions
parser.add_argument('-d','--database',dest='database',
        help="mysql database name or path to sqlite-db [./agres.db]",
        default="./agres.db")
parser.add_argument('-H','--hostname',dest='hostname', 
        help="mysql hostname; if not provided sqlite-db will be used",
        required=False)
parser.add_argument('-u','--user',dest='dbuser', help="mysql database username",
        required=False)
parser.add_argument('-p','--password',dest='mariadbpassword',
        help="mysql password", required=False)
# input specifications
parser.add_argument('-i', '--input_path', dest='input_path',
        help="path to resfinder output dir or tabular input file", required=True)
parser.add_argument('--method', dest='method',
        help="which tools output should be imported", required=True)
parser.add_argument('--mode', dest='mode', 
        help='define which input type was used (fasta/fastq) ' \
                + '- only in use for resfinder [fasta]', default='fasta')

# optional reference to external id or provide a sample_name or both
parser.add_argument('--external_id', dest='external_id', 
        help="optional external id for other db", required=False,
        metavar="INT")
parser.add_argument('--sample_name', dest='sample_name', 
        help="optional sample name to be displayed", required=False,
        metavar="NAME")
parser.add_argument('--assembly', dest='assembly',
        help="imports contig in Fasta-format (dna) to allow visualization",
        required=False, metavar="FASTA")


def main():

    args = parser.parse_args()

    # establish connection to database (mysql or sqlite)
    if args.hostname:
        engine = create_engine("mysql+mysqlconnector://%s:%s@%s:3306/%s" %
                        (args.dbuser, args.mariadbpassword, args.hostname, 
                            args.database))
    else:
        engine = create_engine("sqlite:///%s" % args.database)
    session = Session(engine)
    Base.prepare(engine)

    read_kwargs = {}
    insert_kwargs = {}

    if not args.assembly:
        read_kwargs["extract_coordinates"] = False
    else:
        insert_kwargs["assembly_path"] = args.assembly  

    # interaction with data in database: fetch sample, read results to df and
    # write to database
    sample_args = {}
    if args.sample_name:
        sample_args["name"] = args.sample_name
    if args.external_id:
        sample_args["external_id"] = args.external_id
    if len(sample_args) == 0:
        # create sample without sample name or external id
        associated_sample = get_or_create(session, Sample,
                name=args.sample_name, external_id=args.external_id)
    else:
        associated_sample = get_or_create(session, Sample, **sample_args)
    results_df = read_result(args.input_path, args.method, **read_kwargs)
    insert_into_db(results_df, args.method, associated_sample, session, 
            **insert_kwargs)

    if (args.method == "resfinder" 
            and os.path.exists(f"{args.input_path}/PointFinder_results.txt")):
        pointfinder_df = read_result(args.input_path, "pointfinder", 
                **read_kwargs)
        pointfinder_df["input_type"] = args.mode
        insert_into_db(pointfinder_df, "pointfinder", associated_sample,
                session, **insert_kwargs)

    session.commit()
    session.close()


if  __name__ == "__main__":
    main()
