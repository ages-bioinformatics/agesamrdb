import numpy as np
from Bio import SeqIO

from .util import get_or_create
from .models import ResfinderSequence, Contig, ResfinderResult, PointfinderResult
from sqlalchemy.exc import NoResultFound


def insert_into_resfinder_results(df, associated_sample, session):
    """
    write results to database:
    creates sample and contigs on the fly or retrieves existing
    sequence link is established via accession -> multiples are handled by crc32_hash (of sequence itself)
    if sample feature is not used, only a single sample with no name is ever created
    """
    df = df.replace([np.nan], [None])
    added_results = []
    for i, row in df.iterrows():
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


def insert_into_pointfinder_results(df, associated_sample, session):
    """
    Add Pointfinderresults to database, linked to associated_sample only (no contig possible)
    """
    for i, row in df.iterrows():
        session.add(PointfinderResult(**row, sample_associated=associated_sample))
    session.commit()


def insert_generic_contig_results(df, associated_sample, session, model, to_db_columns, contig_name_col="seqID"):
    """
    write results to database:
    creates sample and contigs on the fly or retrieves existing
    sequence link is established via accession -> multiples are handled by crc32_hash (of sequence itself)
    if sample feature is not used, only a single sample with no name is ever created
    """
    df = df.replace([np.nan], [None])
    added_results = []
    for contig_name, sub_df in df.groupby(contig_name_col):
        associated_contig = session.query(Contig).filter_by(sample_associated=associated_sample, name=contig_name).first()
        if not associated_contig:
            continue

        for i, row in sub_df.iterrows():
            row = row[to_db_columns]
            added_results.append(model(contig_associated=associated_contig, **row))

    session.add_all(added_results)
    session.commit()