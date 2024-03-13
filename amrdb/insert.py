import numpy as np
from Bio import SeqIO

from .util import get_or_create
from .models import ResfinderSequence, Contig, ResfinderResult, \
        PointfinderResult, SpeciesfinderResult, Phenotype, \
        AmrfinderSequence, AmrfinderPointResult, AmrfinderResult
from sqlalchemy.exc import NoResultFound


def _get_linked_sequence_and_contig(row, associated_sample, session, method, contig_kwargs=None):
    if method == "resfinder":
        sequence_model = ResfinderSequence
    elif method == "amrfinder":
        sequence_model = AmrfinderSequence
    if contig_kwargs:
        associated_contig = get_or_create(session, Contig, 
                sample_associated=associated_sample, **contig_kwargs)
    else:
        associated_contig = None
    sequence = session.query(sequence_model).filter_by(accession=row["accession"])
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
        raise NoResultFound("ERROR: missing accession in database: "
                            + f"{row['accession']}") # update database?

    return associated_contig, sequence


def insert_into_resfinder_results(df, associated_sample, session, **kwargs):
    """
    write results to database:
    creates sample and contigs on the fly or retrieves existing
    sequence link is established via accession -> multiples are handled by 
    crc32_hash (of sequence itself) if sample feature is not used, only a
    single sample with no name is ever created
    """
    df = df.replace([np.nan], [None])
    added_results = []
    for i, row in df.iterrows():
        contig_kwargs = {}
        if row.get("contig_name"):
            contig_kwargs["name"] = row["contig_name"]
            contig_kwargs["length"] = row["contig_len"]
        else:
            row["orientation"] = None

        associated_contig, sequence = _get_linked_sequence_and_contig(row,
                associated_sample, session, "resfinder", contig_kwargs)

        row = row[["identity","coverage","ref_pos_start","ref_pos_end",
            "qc_issues","orientation"]]
        added_results.append(ResfinderResult(stored_sequence=sequence,
            contig_associated=associated_contig,
            sample_associated=associated_sample, **row, **kwargs))
    session.add_all(added_results)


def insert_into_amrfinder_results(df, associated_sample, session, **kwargs):
    df = df.replace([np.nan], [None])
    added_results = []
    df_points = df[df["method"].str.startswith("POINT")].copy()
    df = df[~df_points.index].copy()
    for i, row in df.iterrows():
        contig_kwargs = {"name": row["contig_name"]}
        associated_contig, sequence = _get_linked_sequence_and_contig(row,
                associated_sample, session, "amrfinder", contig_kwargs)
        row = row[["identity","coverage","ref_pos_start","ref_pos_end",
                "qc_issues","orientation","method"]]

        added_results.append(AmrfinderResult(stored_sequence=sequence,
            contig_associated=associated_contig,
            sample_associated=associated_sample, **row, **kwargs))

    for i, row in df_points.iterrows():
        contig_kwargs = {"name": row["contig_name"]}
        associated_contig, sequence = _get_linked_sequence_and_contig(row,
                associated_sample, session, "amrfinder", contig_kwargs)
        row = row[["identity","coverage","ref_pos_start","ref_pos_end",
                "qc_issues","orientation","method"]]
        phenotypes = [get_or_create(session, Phenotype, phenotype=p) \
                for p in row["Phenotype"].str.split("/")]
        added_results.append(AmrfinderPointResult(phenotypes=phenotypes,
            contig_associated=associated_contig,
            sample_associated=associated_sample, **row, **kwargs))
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


def add_new_sequences(df, session, tool_model):
    """
    in case the sequence is very similar to known genes and was detected
    but is not 100% identical, we add the entire sequence to be able
    to compare STs of sequences later on and save us from recalculation
    when db updates provide additional STs.
    criteria: no qc issues, identity > 95 < 100, coverage > 95 != 100
    sequence is added with closest related accession and adopts
    these phenotypes (both can be replaced during updates)
    accession is added because in update scenarios, the phenotype should
    be kept also for our added sequences (phenotype association made over acn)
    tool_model: may be ResfinderSequence or AmrfinderSequence
    """
    not_identical = (df["identity"] < 100) | (df["coverage"] != float(100))
    # qc-criteria: no frameshift, start & stopcodon present
    no_issues = (df["qc_issues"].isna())
    above_threshold = (df["identity"] >= 95) & (df["coverage"] > 95)

    for i, row in df[not_identical & no_issues & above_threshold].iterrows():
        entry = session.query(tool_model).filter_by(crc32_hash=row["crc32_hash"])
        if entry.count() > 1: # deal with collisions
            entry = entry.filter_by(sequence=row["sequence"]).first()
        else:
            entry = entry.first()
        if not entry:
            #derive phenotypes from best hit accession - important: Display Warning in UI!
            phenotype_list = session.query(tool_model).filter_by(
                    accession=row["accession"]).first().phenotypes
            session.add(tool_model(
                name=(row["Resistance gene"] + "_AGES_"+row["crc32_hash"]),
                accession=row["accession"], crc32_hash=row["crc32_hash"], 
                sequence=row["sequence"], internal_numbering="AGES_"+row["crc32_hash"],
                phenotypes=phenotype_list))

    session.commit()


def insert_into_pointfinder_results(df, associated_sample, session, **kwargs):
    """
    Add Pointfinderresults to database, linked to associated_sample only
    - there's no contig information for each mutation
    """
    for mutation, sub_df in df.groupby("mutation", as_index=False):
        phenotypes = [get_or_create(session, Phenotype,
            phenotype=p) for p in sub_df["phenotype"].values]
        for i, row in sub_df.iterrows():
            del row["phenotype"]
            session.add(PointfinderResult(**row, sample_associated=associated_sample,
                phenotypes=phenotypes, **kwargs))
    session.commit()
    

def insert_into_speciesfinder_results(df, associated_sample, session, **kwargs):
    """
    Add Speciesfinderresults to database
    """
    for i, row in df.iterrows():
        session.add(SpeciesfinderResult(**row, sample_associated=associated_sample,
            **kwargs))
    session.commit()


def insert_generic_contig_results(df, associated_sample, session, model,
        to_db_columns, contig_name_col="seqID", create_contig=False, **kwargs):
    """
    write results to database, generic function,
    takes a model class as parameter and to_db_columns. model needs to have an
    associated_contig; contigs are not created, entries refering to 
    non-existing contigs will be dropped
    """
    df = df.replace([np.nan], [None])
    for contig_name, sub_df in df.groupby(contig_name_col):
        associated_contig = session.query(Contig).filter_by(
                sample_associated=associated_sample, name=contig_name).first()
        if not associated_contig:
            if create_contig:
                associated_contig = get_or_create(session, Contig, name=contig_name, 
                    sample_associated=associated_sample)
            else:
                continue

        for i, row in sub_df.iterrows():
            row = row[to_db_columns]
            session.add(model(contig_associated=associated_contig, **row, **kwargs))

    session.commit()
