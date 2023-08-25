from Bio import SeqIO
import pandas as pd

from .util import calc_sequence_hash, gene_quality_control

def read_isescan_results(input_file: str) -> pd.DataFrame:
    """
    parses tabular result and return df
    """
    df = pd.read_csv(input_file, sep="\t")
    column_mapping = {"isBegin":"is_start_pos","isEnd":"is_end_pos",
            "ncopy4is":"is_copy_number","start1":"ir_start_pos1",
            "end1":"ir_end_pos1","start2":"ir_start_pos2","end2":"ir_end_pos2",
            "orfBegin":"orf_start_pos","orfEnd":"orf_end_pos",
            "strand":"orientation", "E-value":"e_value",
            }
    df = df.rename(columns=column_mapping)
    df["complete"] = (df["type"] == "c")

    return df


def read_bakta_results(input_file: str) -> pd.DataFrame:
    """
    parses tabular result and return df
    """
    header = ["seqID","product_type","ref_pos_start","ref_pos_end",
            "orientation","locus_tag","gene_name","product","db_xref"]
    df = pd.read_csv(input_file, sep="\t", comment="#", names=header)

    return df


def read_mobtyper_results(input_file: str) -> pd.DataFrame:
    """
    parses tab-separated result and return df
    """
    df = pd.read_csv(input_file, sep="\t")
    column_mapping = {c: c.replace("(s)","") for c in df.columns if "(s)" in c}
    column_mapping["gc"] = "gc_content"
    df["seqID"] = df["sample_id"].str.split(" ", expand=True)[0]
    df = df.rename(columns=column_mapping)

    return df
    

def read_resfinder_results(resfinder_dir: str, 
        extract_coordinates: bool=True) -> pd.DataFrame:
    """
    parses sequence hits and tabular result to extract all needed information
    extract_coordinates: boolean flag indicates wheter coordinates can be
        extracted
    """
    # read sequences to df
    sequences = []
    sequence_file=f"{resfinder_dir}/ResFinder_Hit_in_genome_seq.fsa"
    for i, seqrecord in enumerate(SeqIO.parse(sequence_file, "fasta")):
        comment = gene_quality_control(seqrecord)
        sequences.append((seqrecord.id.strip(","), str(seqrecord.seq), comment,
            calc_sequence_hash(str(seqrecord.seq)), seqrecord.description))

    columns = ["Resistance gene", "sequence", "qc_issues","crc32_hash","desc"]
    df_sequences = pd.DataFrame(sequences, columns=columns)

    pos_regex = "Contig name: (.*?), Position: ([NA0-9]*\.\.[NA0-9]*)"
    df_sequences[["Contig","Position in contig"]] = df_sequences["desc"].str.extract(pos_regex, expand=True)

    # read tabluar resfinder output to df
    df_results = pd.read_csv(f"{resfinder_dir}/ResFinder_results_tab.txt",
            sep="\t")

    # handle case that no resistance gene was found
    if df_results.empty:
        df_results[["ref_pos_end","ref_pos_start","qc_issues"]] = None
        df = df_results
    # merge seq and tabluar resfinder output:
    elif extract_coordinates:
        # on 3 pseudo unique positions and drop duplicates
        # duplicate case: multiple accessions (alleles) matching equally good at same pos in contig
        # drop duplicates on contig and position in contig
        df = df_results.merge(df_sequences, on=["Resistance gene","Contig",
            "Position in contig"]).drop_duplicates(["Contig","Position in contig"])
        df["Contig"] = df["Contig"].str.split(" ", expand=True)[0]
        df[["ref_pos_start","ref_pos_end"]] = df["Position in contig"].str.split("..", regex=False, expand=True)
        df["ref_pos_start"] = df["ref_pos_start"].astype(int)
        df["ref_pos_end"] = df["ref_pos_end"].astype(int)
    else:
        # on resistance gene name only. drop duplicates on sequence hash
        # duplicate case: multiple alleles matching equally good, no posistion
        # extractable -> remove duplicate sequences
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


def read_pointfinder_results(resfinder_dir: str) -> pd.DataFrame:
    """
    parses pointfinder result and return df
    """
    df = pd.read_csv(f"{resfinder_dir}/PointFinder_results.txt", sep="\t")
    df["phenotype"] = df["Resistance"].str.split(",")
    df = df.explode("phenotype")
    df["phenotype"] = df["phenotype"].apply(lambda k: k.strip().title())
    df = df.drop_duplicates(["Mutation","phenotype"])
    df.rename(columns={"Mutation":"mutation","Nucleotide change":"nuc_change"},
            inplace=True)
    df = df[["phenotype","mutation","nuc_change"]]
    return df

