#!/usr/bin/env python

import sys
import os
import glob
import argparse

import numpy as np
import pandas as pd

from Bio import SeqIO

from amrdb.models import Base, ResfinderSequence, Phenotype, ResfinderResult, \
        phenotype_association_table, Sample, Contig, PointfinderResult, \
        ISEScanResult, BaktaResult, MobTyperResult
from amrdb.util import calc_sequence_hash

from sqlalchemy import create_engine, insert, text, inspect 
from sqlalchemy.orm import Session, declarative_base
from sqlalchemy.ext.automap import automap_base

parser = argparse.ArgumentParser(description="Create and/or update Resistance" \
        + "-Gene-Database for ResFinder")
parser.add_argument('-d','--database',dest='database',
        help="mysql database name or path to sqlite-db [./agres.db]",
        default="./agres.db")
parser.add_argument('-H','--hostname',dest='hostname',
        help="mysql hostname if not provided sqlite-db will be created",
        required=False)
parser.add_argument('-u','--user',dest='dbuser', help="mysql database username",
        required=False)
parser.add_argument('-p','--password',dest='mariadbpassword',
        help="mysql password", required=False)
parser.add_argument('--resfinder_db', dest='resfinder_db_dir',
        help="path to resfinder github checkout", required=False,
        default="resfinder_db")


def read_resfinder_databases(resfinder_db_dir):
    """
    parses all fasta-files and stores in record dict of records (list of dicts)
    """
    entries = []
    for db_file in glob.glob(resfinder_db_dir+"/*.fsa"):
        if db_file.endswith("all.fsa"):
            continue
        table_name = os.path.basename(db_file).replace(".fsa","").replace("-","_")
        for seqrecord in SeqIO.parse(db_file, "fasta"):
            entry = {"name": seqrecord.id, "sequence":str(seqrecord.seq), 
                    "crc32_hash":calc_sequence_hash(str(seqrecord.seq))}
            entries.append(entry)
        
    return entries


def read_phenotypes(resfinder_db_dir):
    name_split_regex = "([^_]*)[-_]([^_]*)_([A-Z]+_[A-Z0-9.]*|[A-Z0-9.]*):*.*$"
    df = pd.read_csv(f"{resfinder_db_dir}/phenotypes.txt", sep="\t")
    df["Phenotype"] = df["Phenotype"].str.split(", ")
    df = df.explode("Phenotype")
    df[["short_name","subseq_numbering","sequence_identifier"]] = \
            df["Gene_accession no."].str.extract(name_split_regex, expand=True)
    phenotypes = df.drop_duplicates(["Phenotype","Class"])[["Phenotype","Class"]]
    phenotypes = phenotypes[~phenotypes["Class"].str.contains(",")]
    phenotypes = phenotypes.reset_index().drop("index",axis=1)
    phenotypes = phenotypes.merge(df[["Phenotype","sequence_identifier"]], 
            on="Phenotype")
    return phenotypes


def identify_columns(records):
    """
    uses records list of dicts to create dataframe and splits name using regex
        into short_name, main_numbering, subseq_numbering and accession
    subseq_numbering: internally used as allels of z.B. sul1 which are never
        shown to user (e.g. 1)
    accession: important unique accession of external reference sequence 
    main_numbering: number of resistance gene variant that is displayed to
        user (e.g. 52B)
    short_name: basic gene name (e.g. blaTEM)
    """    
    name_split_regex = "([^_]*)[-_]([^_]*)_([A-Z]+_[A-Z0-9.]*|[A-Z0-9.]*):*.*$"
    df = pd.DataFrame.from_records(records)
    df[["short_name","subseq_numbering","accession"]] = \
            df["name"].str.extract(name_split_regex, expand=True)

    try:
        df[["short_name", "main_numbering"]] = \
                df["short_name"].str.extract("([^\d^(]*)-*([()0-9A-Za-z'-]*)", expand=True)
    except ValueError:
        df["main_numbering"] = ""

    return df


def update_existing_sequences(df_new, session):
    """
    Updates all sequence-names, numbering and accession in database:
    not sequence, crc32_hash or internal identifier unless newly added
    bc. already found variants/alleles are added with internal identifier
    and closest accession (single acn. may have multiple sequences associated)
    this saves some recalculation for newly added variants which we already
    detected in our assemblies (just adding name)
    """
    df_new.replace([np.nan], [None], inplace=True)
    for i, row in df_new.iterrows():
        entry = session.query(ResfinderSequence).filter_by(
                crc32_hash=row["crc32_hash"])
        if entry.count() > 1: # deal with collisions
            entry = entry.filter_by(sequence=row["sequence"]).first()
        else:
            entry = entry.first()
        if entry:
            print(f"updating {entry.name} with {row}")
            for key, value in row.items():
                setattr(entry, key, value)
        else:
            #TODO: handle that samples should be reanalyzed
            # potentially to be checked manually wheter required or not
            # i.e. only in case a completely new type of sequence 
            # not variant/allel is added
            print(f"adding new {row}")
            session.add(ResfinderSequence(**row))

    session.commit()


def write_phenotypes(phenotypes_df, session):
    """
    create entries for phenotypes adding links to sequence by accession
    by recreating from scratch we update the phenotype associated with sequences
    e.g. previously added internal sequence receives new phenotype association
    if they are published in database with official name
    """
    for groupname, sub_df in phenotypes_df.groupby(["Phenotype","Class"]):
        linked_sequences = session.query(ResfinderSequence).filter(
                ResfinderSequence.accession.in_(
                    sub_df["sequence_identifier"].to_list())).all()
        session.add(Phenotype(phenotype=groupname[0], class_name=groupname[1],
            sequences=linked_sequences))
    session.commit()


def main():

    args = parser.parse_args()

    # read database files:
    records = read_resfinder_databases(args.resfinder_db_dir)
    sequences_df = identify_columns(records)
    phenotypes_df = read_phenotypes(args.resfinder_db_dir)

    # establish connection to database
    if args.hostname:
        engine = create_engine("mysql+mysqlconnector://%s:%s@%s:3306/%s" %
                        (args.dbuser, args.mariadbpassword, args.hostname,
                            args.database))
    else:
        engine = create_engine("sqlite:///%s" % args.database)

    session = Session(engine)

    # get all tables from database
    insp = inspect(engine)
    if not insp.has_table(ResfinderSequence.__table__.name):
        # create tables if not existing yet, only first initialization
        # assume others also not existing (i.e when run for first time)
        ResfinderSequence.__table__.create(engine) 
        Sample.__table__.create(engine)
        Contig.__table__.create(engine)
        ResfinderResult.__table__.create(engine)
        PointfinderResult.__table__.create(engine)
        ISEScanResult.__table__.create(engine)
        BaktaResult.__table__.create(engine)
        MobTyperResult.__table__.create(engine)

    # always drop phenotype table/association table and create newly
    if insp.has_table(phenotype_association_table.name):
        phenotype_association_table.drop(engine)
    if insp.has_table(Phenotype.__table__.name):
        Phenotype.__table__.drop(engine)
    Phenotype.__table__.create(engine)
    phenotype_association_table.create(engine)
    session.commit()    

    # connect the defined classes to the data in db:
    Base.prepare(engine)

    # write database entries:
    update_existing_sequences(sequences_df, session)
    write_phenotypes(phenotypes_df, session)

    session.commit()
    session.close()


if  __name__ == "__main__":
    main()
