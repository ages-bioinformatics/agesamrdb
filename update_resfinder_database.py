#!/usr/bin/env python

import sys
import os
import glob
import argparse

import numpy as np
import pandas as pd

from Bio import SeqIO

import agesamrdb
from agesamrdb.models import Base, ResfinderSequence, Phenotype, ResfinderResult, \
        Sample, Contig, PointfinderResult, ISEScanResult, BaktaResult, \
        MobTyperResult, InVitroResult, PlasmidfinderResult, PhispyResults, \
        SpeciesfinderResult, pointfinder_phenotype_association_table, ToolVersion, \
        AmrfinderSequence, AmrfinderResult, AmrfinderPointResult, \
        amrfinder_point_phenotype_association_table, \
        amrfinder_phenotype_association_table, resfinder_phenotype_association_table

from agesamrdb.util import calc_sequence_hash, get_or_create

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
        help="path to resfinder_db local git checkout", required=False,
        default="resfinder_db")
#parser.add_argument('--pointfinder_db', dest='pointfinder_db_dir',
#        help="path to pointfinder_db local git checkout (currently not functional)",
#        required=False, default='pointfinder_db')
parser.add_argument('--amrfinder_db', dest='amrfinder_db_dir',
        help="path to amrfinder_db previously downloaded",
        required=False, default='/db/amrfinder/latest')


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
        
    return pd.DataFrame.from_records(entries)


def read_amrfinder_database(amrfinder_db_dir):
    # TODO read table for families too - currently there are 2000+ sequences
    # which have no phenotype associated (because part of families)
    entries = []
    db_file = os.path.join(amrfinder_db_dir,"AMRProt")
    for seqrecord in SeqIO.parse(db_file, "fasta"):
        protein_gi, accession, _, _, name, short_name, activity_type, core_status,\
                phenotype, class_name, long_name = seqrecord.description.split("|")
        entry = {"name": name, "sequence":str(seqrecord.seq),
                "crc32_hash":calc_sequence_hash(str(seqrecord.seq)),
                "short_name": short_name, "activity_type": activity_type,
                "Phenotype": phenotype, "Class": class_name, "long_name": long_name,
                "is_core": True if int(core_status) > 1 else False,
                "accession": accession,
                }
        entries.append(entry)
    return pd.DataFrame.from_records(entries)


def read_phenotypes(resfinder_db_dir):
    name_split_regex = "([^_]*)[-_]([^_]*)_([A-Z]+_[A-Z0-9.]*|[A-Z0-9.]*):*.*$"
    df = pd.read_csv(f"{resfinder_db_dir}/phenotypes.txt", sep="\t")
    df["Phenotype"] = df["Phenotype"].str.split(",")
    df = df.explode("Phenotype")
    df["Phenotype"] = df["Phenotype"].str.title().str.replace("_", " ").str.strip()
    df[["short_name","subseq_numbering","sequence_identifier"]] = \
            df["Gene_accession no."].str.extract(name_split_regex, expand=True)
    return df


def read_pointfinder_phenotypes(pointfinder_db_dir):
    """
    read all pointmutations from database and establish connection between
    phenotypes and mutations. problem: pointfinder reports mutation differently
    and generates the output name in very complex decision tree.
    #TODO write a function that imports all possible pointfinder results before
    #current behaviour: create database entries on the fly when result is imported
    """
    #TODO finish this function
    #
    search_pattern = f"pointfinder_db_dir/*/resistens-overview.txt"
    for phenotype_table_file in glob.glob(search_pattern):
        organism = os.path.basename(os.path.dirname(phenotype_table_file))
        df = pandas.read_csv(phenotype_table_file, sep="\t")
        #


def identify_columns(df):
    """
    uses dataframe and splits name using regex
        into short_name, main_numbering, subseq_numbering and accession
    subseq_numbering: internally used as allels of z.B. sul1 which are never
        shown to user (e.g. 1)
    accession: important unique accession of external reference sequence 
    main_numbering: number of resistance gene variant that is displayed to
        user (e.g. 52B)
    short_name: basic gene name (e.g. blaTEM)
    """    
    name_split_regex = "([^_]*)[-_]([^_]*)_([A-Z]+_[A-Z0-9.]*|[A-Z0-9.]*):*.*$"
    df[["short_name","subseq_numbering","accession"]] = \
            df["name"].str.extract(name_split_regex, expand=True)

    try:
        df[["short_name", "main_numbering"]] = \
                df["short_name"].str.extract("([^\d^(]*)-*([()0-9A-Za-z'-]*)", expand=True)
    except ValueError:
        df["main_numbering"] = ""

    return df


def update_existing_sequences(df_new, session, tool_model):
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
        entry = session.query(tool_model).filter_by(
                crc32_hash=row["crc32_hash"])
        if entry.count() > 1: # deal with collisions
            entry = entry.filter_by(sequence=row["sequence"]).first()
        else:
            entry = entry.first()
        if entry:
            print(f"updating {entry} with {row}")
            for key, value in row.items():
                setattr(entry, key, value)
        else:
            #TODO: handle that samples should be reanalyzed
            # potentially to be checked manually wheter required or not
            # i.e. only in case a completely new type of sequence 
            # not variant/allel is added
            print(f"adding new {row}")
            session.add(tool_model(**row))

    session.commit()


def write_phenotypes(phenotypes_df, session, tool_model):
    """
    create entries for phenotypes adding links to sequence by accession
    by recreating from scratch we update the phenotype associated with sequences
    e.g. previously added internal sequence receives new phenotype association
    if they are published in database with official name
    """
    for groupname, sub_df in phenotypes_df.groupby("Phenotype"):
        phenotype_name = groupname.title().replace("_"," ")
        linked_sequences = session.query(tool_model).filter(
                tool_model.accession.in_(
                    sub_df["sequence_identifier"].to_list())).all()
        phenotype = get_or_create(session, Phenotype, phenotype=phenotype_name)
        for l in linked_sequences:
            if tool_model.__name__.startswith("Resfinder"):
                phenotype.resfinder_sequences.append(l)
            elif tool_model.__name__.startswith("Amrfinder"):
                phenotype.amrfinder_sequences.append(l)
        session.add(phenotype)
    session.commit()


def initialize_phenotype_classes(session):
    """
    initialize phenotype - class table, manually curated to avoid wrong class
    /substance associations, which are not readable automatically without errors
    from resfinder or amrfinder.
    TODO: might be better curated though
    """
    pkg_dir = os.path.dirname(os.path.abspath(agesamrdb.__file__))
    phenotype_table_csv = os.path.join(pkg_dir,"data","phenotypes_classnames.tsv")
    phenotype_df = pd.read_csv(phenotype_table_csv, sep="\t")
    for i, row in phenotype_df.iterrows():
        phenotype = get_or_create(session, Phenotype, **row)
        session.add(phenotype)
    session.commit()


def main():

    args = parser.parse_args()

    # read database files:
    records = read_resfinder_databases(args.resfinder_db_dir)
    sequences_df = identify_columns(records)
    phenotypes_df = read_phenotypes(args.resfinder_db_dir)
    amrfinder_df = read_amrfinder_database(args.amrfinder_db_dir)
    amrfinder_phenotype_df = amrfinder_df[["Phenotype","accession"]]\
            .rename(columns={"accession":"sequence_identifier"})
    amrfinder_phenotype_df["Phenotype"] = amrfinder_phenotype_df["Phenotype"].str.split("/")
    amrfinder_phenotype_df = amrfinder_phenotype_df.explode("Phenotype")
    amrfinder_phenotype_df["Phenotype"] = amrfinder_phenotype_df["Phenotype"].str.title()
    amrfinder_df = amrfinder_df.drop(["Phenotype","Class"], axis=1)

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
    install = False
    if not insp.has_table(ResfinderSequence.__table__.name):
        print("create new database")
        # create tables if not existing yet, only first initialization
        # assume others also not existing (i.e when run for first time)
        install = True
        ResfinderSequence.__table__.create(engine) 
        AmrfinderSequence.__table__.create(engine)
        Sample.__table__.create(engine)
        Contig.__table__.create(engine)
        ToolVersion.__table__.create(engine)
        ResfinderResult.__table__.create(engine)
        AmrfinderResult.__table__.create(engine)
        AmrfinderPointResult.__table__.create(engine)
        PointfinderResult.__table__.create(engine)
        #ISEScanResult.__table__.create(engine)
        BaktaResult.__table__.create(engine)
        MobTyperResult.__table__.create(engine)
        Phenotype.__table__.create(engine)
        InVitroResult.__table__.create(engine)
        PlasmidfinderResult.__table__.create(engine)
        ISEScanResult.__table__.create(engine)
        PhispyResults.__table__.create(engine)
        SpeciesfinderResult.__table__.create(engine)
        pointfinder_phenotype_association_table.create(engine)
        amrfinder_point_phenotype_association_table.create(engine)
        amrfinder_phenotype_association_table.create(engine)
    
    # always drop phenotype association table and create newly
    if insp.has_table(resfinder_phenotype_association_table.name):
        resfinder_phenotype_association_table.drop(engine)
    if insp.has_table(amrfinder_phenotype_association_table.name):
        amrfinder_phenotype_association_table.drop(engine)

    amrfinder_phenotype_association_table.create(engine)
    resfinder_phenotype_association_table.create(engine)
    session.commit()

    # connect the defined classes to the data in db:
    Base.prepare(engine)

    if install:
        initialize_phenotype_classes(session)

    # write database entries:
    print("start filling resfinder sequences to database")
    update_existing_sequences(sequences_df, session, ResfinderSequence)
    write_phenotypes(phenotypes_df, session, ResfinderSequence)
    print("start filling amrfinder sequences to database")
    update_existing_sequences(amrfinder_df, session, AmrfinderSequence)
    write_phenotypes(amrfinder_phenotype_df, session, AmrfinderSequence)

    session.commit()
    session.close()


if  __name__ == "__main__":
    main()
