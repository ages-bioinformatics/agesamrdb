#!/usr/bin/env python

import sys
import os
import glob
import argparse

import numpy as np
import pandas as pd

from Bio import SeqIO

from amrdb.models import Base, ResfinderSequence, Phenotype, ResfinderResult, Sample, Contig
from amrdb.util import calc_sequence_hash, gene_quality_control

from sqlalchemy import create_engine, inspect, select
from sqlalchemy.orm import Session
from sqlalchemy.exc import NoResultFound

parser = argparse.ArgumentParser(description="Add single Refinder result to amrdb")
parser.add_argument('-d','--database',dest='database', help="mysql database name", required=True)
parser.add_argument('-H','--hostname',dest='hostname', help="mysql hostname", required=True)
parser.add_argument('-u','--user',dest='dbuser', help="mysql database username", required=True)
parser.add_argument('-p','--password',dest='mariadbpassword', help="mysql password", required=True)
parser.add_argument('-i', '--resfinder_output', dest='resfinder_dir', help="path to resfinder output dir", required=True)
parser.add_argument('--external_id', dest='external_id', help="optional external id for other db", required=False)
parser.add_argument('--sample_name', dest='sample_name', help="optional sample name to be displayed", required=False)
parser.add_argument('--mode', dest='mode', help="mode wether fastq or fasta input", required=False, default="fasta")
parser.add_argument('--assembly', dest='assembly', help="allows contig visualization", required=False)



def read_resfinder_results(resfinder_dir, extract_coordinates=True):
    """
    parses sequence hits and tabular result to extract all needed information
    extract_coordinates: boolean flag indicates wheter coordinates can be extracted
    """
    # read sequences to df
    sequences = []
    for i, seqrecord in enumerate(SeqIO.parse(f"{resfinder_dir}/ResFinder_Hit_in_genome_seq.fsa", "fasta")):
        comment = gene_quality_control(seqrecord)
        sequences.append((seqrecord.id.strip(","), str(seqrecord.seq), comment, calc_sequence_hash(str(seqrecord.seq)), seqrecord.description))

    df_sequences = pd.DataFrame(sequences, columns=["Resistance gene", "sequence", "qc_issues","crc32_hash","desc"])
    df_sequences[["Contig","Position in contig"]] = df_sequences["desc"].str.extract("Contig name: (.*?), Position: ([NA0-9]*\.\.[NA0-9]*)", expand=True)

    # read tabluar resfinder output to df
    df_results = pd.read_csv(f"{resfinder_dir}/ResFinder_results_tab.txt", sep="\t")

    # merge seq and tabluar resfinder output:
    if extract_coordinates:
        # on 3 pseudo unique positions and drop duplicates
        # duplicate case: multiple accessions (alleles) matching equally good at same pos in contig
        # drop duplicates on contig and position in contig
        df = df_results.merge(df_sequences, on=["Resistance gene","Contig","Position in contig"]).drop_duplicates(["Contig","Position in contig"])
        df["Contig"] = df["Contig"].str.split(" ", expand=True)[0]
        df[["ref_pos_start","ref_pos_end"]] = df["Position in contig"].str.split("..", regex=False, expand=True)
        df["ref_pos_start"] = df["ref_pos_start"].astype(int)
        df["ref_pos_end"] = df["ref_pos_end"].astype(int)
    else:
        # on resistance gene name only. drop duplicates on sequence hash
        # duplicate case: multiple alleles matching equally good, no posistion extractable
        # remove duplicate sequences
        df = df_results.merge(df_sequences, on="Resistance gene").drop_duplicates("crc32_hash")
        df["ref_pos_end"] = None
        df["ref_pos_start"] = None


    df.rename(columns={"Accession no.":"accession",
        "Identity":"identity",
        "Coverage":"coverage",
        "Contig":"contig_name",
        },
        inplace=True,
    )

    return df


def add_new_sequences(df, session):
    """
    in case the sequence is very similar to known genes and was detected
    but is not 100% identical, we add the entire sequence to be able
    to compare STs of sequences later on and save us from recalculation
    when db updates provide additional STs.
    criteria: no qc issues, identity > 95 < 100, coverage > 95 != 100
    sequence is added with closest related accession and adopts
    these phenotypes (both can be replaced during updates)
    """
    not_identical = (df["identity"] < 100) | (df["coverage"] != float(100))
    no_issues = (df["qc_issues"].isna()) # qc-criteria: no frameshift, start & stopcodon present
    above_threshold = (df["identity"] >= 95) & (df["coverage"] > 95) # e.g. in-frame deletions

    for i, row in df[not_identical & no_issues & above_threshold].iterrows():
        entry = session.query(ResfinderSequence).filter_by(crc32_hash=row["crc32_hash"])
        if entry.count() > 1: # deal with collisions
            entry = entry.filter_by(sequence=row["sequence"]).first()
        else:
            entry = entry.first()
        if not entry:
            #derive phenotypes from best hit accession - important: Display Warning in UI!
            phenotype_list = session.query(ResfinderSequence).filter_by(accession=row["accession"]).first().phenotypes
            session.add(ResfinderSequence(name=(row["Resistance gene"]+"_AGES_"+row["crc32_hash"]),
                sequence=row["sequence"], accession=row["accession"], crc32_hash=row["crc32_hash"],
                internal_numbering="AGES_"+row["crc32_hash"], phenotypes=phenotype_list))

    session.commit()


def add_contig_info(df, assembly_file):
    """
    if applicable, parses assembly and reads orientation of detected gene
    together with total length of contig (for displaying purposes)
    """
    contig_names = df["contig_name"].unique()
    contigs = {}
    contig_lens = {}
    for record in SeqIO.parse(assembly_file, "fasta"):
        if record.id in contig_names:
            contigs[record.id] = record
            contig_lens[record.id] = len(record)

    df["contig_len"] = df["contig_name"].map(contig_lens)

    ori = {}
    for contig_name, sub in df.groupby("contig_name"):
        for i, row in sub.iterrows():
            seq = contigs[contig_name][row["ref_pos_start"]-1:row["ref_pos_end"]].seq.replace("-","")
            sequence = row["sequence"].replace("-","")
            if str(seq) == sequence:
                ori[i] = "+"
            elif str(seq.reverse_complement()) == sequence:
                ori[i] = "-"
            else:
                # should never occur
                ori[i] = np.nan
                print([i, "mismatch", seq, seq.reverse_complement(), sequence])

    df["orientation"] = ori
    return df


def get_or_create(session, model, **kwargs):
    """
    generic function similar to what's known from django ORM
    may throw an exeption when kwargs would match more than one sample
    if only get_or_create used to create objects, this will never happen
    """
    instance = session.query(model).filter_by(**kwargs).one_or_none()
    if instance:
        return instance
    else:
        instance = model(**kwargs)
        session.add(instance)
        session.commit()
        return instance


def insert_resfinder_results(df, session):
    """
    write results to database:
    creates sample and contigs on the fly or retrieves existing
    sequence link is established via accession -> multiples are handled by crc32_hash (of sequence itself)
    if sample feature is not used, only a single sample with no name is ever created
    """
    df = df.replace([np.nan], [None])
    added_results = []
    for i, row in df.iterrows():
        associated_sample = get_or_create(session, Sample, name=row["sample_name"], external_id=row["external_id"])
        if row.get("contig_name"):
            associated_contig = get_or_create(session, Contig, length=row["contig_len"], name=row["contig_name"], sample_associated=associated_sample)
        else:
            associated_contig = None
            row["orientation"] = None
        sequence = session.query(ResfinderSequence).filter_by(accession=row["accession"])
        if sequence.count() > 1:
            exact_matches = sequence.filter_by(crc32_hash=row["crc32_hash"])
            if exact_matches.count() != 1:
                # no perfect match -> use first with accession (original entry)
                sequence = sequence.first()
            else:
                sequence = exact_matches.first()
        else:
            sequence = sequence.first()
        if not sequence:
            raise NoResultFound(f"ERROR: missing accession in database: {row['accession']}") # update database?
        row = row[["identity","coverage","ref_pos_start","ref_pos_end","qc_issues","orientation"]]
        added_results.append(ResfinderResult(stored_sequence=sequence, contig_associated=associated_contig, sample_associated=associated_sample, **row))
    session.add_all(added_results)
    session.commit() 


def main():

    args = parser.parse_args()
    if args.mariadbpassword and args.dbuser:
        mariadbpassword = args.mariadbpassword
        dbuser = args.dbuser

    extract_coordinates = False
    if args.mode == "fasta":
        extract_coordinates = True

    results_df = read_resfinder_results(args.resfinder_dir, extract_coordinates)
    
    results_df["external_id"] = args.external_id
    results_df["sample_name"] = args.sample_name
    
    if args.assembly:
        results_df = add_contig_info(results_df, args.assembly)

    # establish connection to database
    engine = create_engine("mysql+mysqlconnector://%s:%s@%s:3306/%s" %
                          (dbuser, mariadbpassword, args.hostname, args.database))
    session = Session(engine)

    # get all tables from database
    # maybe not needed after all, as they need to exist always before Base.prepare
    # --> are created at update script
    #insp = inspect(engine)
    #if not insp.has_table(ResfinderResult.__table__.name):
    #    # create table if not existing yet, only first initialization
    #    ResfinderResult.__table__.create(engine) 
    #    session.commit()    

    Base.prepare(engine)
    if args.mode == "fasta":
        # currently only allow new gene entries from assemblies to avoid
        # incomplete gap-containing constructs from reads or
        # spamming of multiple similar variants caused by sequencing errors
        add_new_sequences(results_df, session)

    insert_resfinder_results(results_df, session)

    session.commit()
    session.close()


if  __name__ == "__main__":
    main()
