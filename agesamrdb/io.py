from Bio import SeqIO
import os
import pandas as pd
import json

from .util import calc_sequence_hash, gene_quality_control, \
    apply_offset_to_partial_contigs

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
    df["seqID"] = df["seqID"].astype(str)
    df = apply_offset_to_partial_contigs(df)
    return df


def read_bakta_results(input_file: str) -> pd.DataFrame:
    """
    parses tabular result and return df
    """
    header = ["seqID","product_type","ref_pos_start","ref_pos_end",
            "orientation","locus_tag","gene_name","product","db_xref"]
    df = pd.read_csv(input_file, sep="\t", comment="#", names=header)
    df["seqID"] = df["seqID"].astype(str)
    df = apply_offset_to_partial_contigs(df)
    return df


def read_mobtyper_results(input_file: str) -> pd.DataFrame:
    """
    parses tab-separated result and return df
    """
    df = pd.read_csv(input_file, sep="\t")
    column_mapping = {c: c.replace("(s)","") for c in df.columns if "(s)" in c}
    column_mapping["gc"] = "gc_content"
    df["seqID"] = df["seqID"].astype(str)
    df["seqID"] = df["sample_id"].str.split(" ", expand=True)[0]
    df = df.rename(columns=column_mapping)

    return df
    

def parse_fasta_hits(sequence_file: str, translate: bool=False) -> pd.DataFrame:
    """
    parses sequence fasta originating from hits (and might translate)
    """
    sequences = []
    for i, seqrecord in enumerate(SeqIO.parse(sequence_file, "fasta")):
        # amrfinder does not include stop-codon
        comment = gene_quality_control(seqrecord, ignore_missing_stop=translate)
        if not comment and translate:
            seq = str(seqrecord.seq.translate())
        else:
            seq = str(seqrecord.seq)
        sequences.append((seqrecord.id.strip(","), seq, comment,
            calc_sequence_hash(seq), seqrecord.description))

    columns = ["Resistance gene", "sequence", "qc_issues","crc32_hash","desc"]
    df_sequences = pd.DataFrame(sequences, columns=columns)

    return df_sequences


def read_resfinder_results(resfinder_dir: str, 
        extract_coordinates: bool=True) -> pd.DataFrame:
    """
    parses sequence hits and tabular result to extract all needed information
    extract_coordinates: boolean flag indicates wheter coordinates can be
        extracted
    """
    # read sequences to df
    df_sequences = parse_fasta_hits(os.path.join(resfinder_dir,"ResFinder_Hit_in_genome_seq.fsa"))

    pos_regex = "Contig name: (.*?), Position: ([NA0-9]*\.\.[NA0-9]*)"
    df_sequences[["Contig","Position in contig"]] = df_sequences["desc"].str.extract(pos_regex, expand=True)

    # read tabluar resfinder output to df
    df_results = pd.read_csv(f"{resfinder_dir}/ResFinder_results_tab.txt",
            sep="\t")
    df_results["Contig"] = df_results["Contig"].astype(str)

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
    df["phenotype"] = df["phenotype"].apply(lambda k: k.strip().replace("_", " ").title())
    df = df.drop_duplicates(["Mutation","phenotype"])
    df.rename(columns={"Mutation":"mutation","Nucleotide change":"nuc_change"},
            inplace=True)
    df = df[["phenotype","mutation","nuc_change"]]
    return df

def read_plasmidfinder_results(input_file: str) -> pd.DataFrame:
    column_mapping = {'Database': 'database_name', 'Plasmid': 'plasmid', 'Identity': 'identity', 
                        'Note': 'note', 'Accession number': 'accession_number', 'Contig': 'contig_name' }
    
    columns = ['database_name', 'plasmid', 'identity', 'contig_name', 'note', 'accession_number', 
        'query_length', 'template_length', 'ref_pos_start', 'ref_pos_end']

    
    df = pd.read_csv(input_file, sep='\t')
    if df.empty:
        df = pd.DataFrame()
        df[columns] = None
    else:
        df[['query_length', 'template_length']] = df['Query / Template length'].str.split(' / ', n=1, expand=True).astype(int)
        df[['ref_pos_start', 'ref_pos_end']] = df['Position in contig'].str.split('\.\.', n=1, expand=True).astype(int)
        df['Contig'] = df['Contig'].astype(str).str.split(" ", expand=True)[0]
        df = df. rename(columns=column_mapping)
        df = df[columns]
    return df
    
def read_phispy_results(input_file: str) -> pd.DataFrame:
    columns = ['prophage_number', 'contig_name', 'start' ,'stop', 'start_attL', 'end_attL', 
                'start_attR', 'end_attR', 'sequence_attL', 'sequence_attR', 'description']
    df = pd.read_csv(input_file, index_col=False, header=None, names=columns, sep='\t')
    df["contig_name"] = df["contig_name"].astype(str)
    print(df)
    df = apply_offset_to_partial_contigs(df, "contig_name")
    return df 
    

def read_speciesfinder_results(input_file: str) -> pd.DataFrame:
    """
    parses CGE Speciesfinder result (json str) and return df
    """
    column_mapping = { 'Template':'template',
                       'Species':'species',
                       'Match':'match_id',
                       'Database':'database_name',
                       'Confidence of result':'confidence_of_result',
                     }

    with open(input_file) as json_file:
        json_data = json.load(json_file)
    
    df = pd.DataFrame.from_records([json_data['speciesfinder']['results']])
    df[['file_format', 'method']] = json_data['speciesfinder']['user_input']['file_format'], \
                                    json_data['speciesfinder']['user_input']['method']
    df = df.rename(columns=column_mapping)
    return df


def read_amrfinder_results(input_dir: str) -> pd.DataFrame:
    """
    important! default output for amrfinder is only tabular file
    input required here is an output directory which contains two
    files: "amrfinder_results.txt" (-o/--output) and
    "amrfinder_nucleotides.fasta" (--nucleotide_output)
    both are assumed to be present in the input directoy (and named accordingly)
    !!!this is not a default output-naming scheme of amrfinder, only our workflow!!!
    """
    # renaming scheme for tabular result (and used columns)
    column_mapping = {
            "Contig id": "contig_name",
            "Start": "ref_pos_start",
            "Stop": "ref_pos_end",
            "Strand": "orientation",
            "Gene symbol": "name",
            "Sequence name": "long_name",
            "Scope": "scope",
            "Method": "method",
            "Subclass": "Phenotype",
            "% Coverage of reference sequence": "coverage",
            "% Identity to reference sequence": "identity",
            "Accession of closest sequence": "accession",
    }
    # parse tabular result
    df = pd.read_csv(os.path.join(input_dir,"amrfinder_results.txt"), sep="\t")
    df = df.rename(columns = column_mapping)
    df = df[[v for v in column_mapping.values()]].copy()
    df["is_core"] = (df["scope"] == "core")
    df["mutation"] = df["name"]
    df["contig_name"] = df["contig_name"].astype(str)

    # parse fasta-result
    df_sequences = parse_fasta_hits(os.path.join(input_dir,"amrfinder_nucleotides.fasta"),
            translate=True)

    # extract position-on-contig and merge to tabular result
    pos_regex = r"^(.*?):(\d+)-(\d+).*([+-]) .*"
    df_sequences[["contig_name","ref_pos_start","ref_pos_end","orientation"]] =\
            df_sequences["desc"].str.extract(pos_regex, expand=True)

    df_sequences["ref_pos_start"] = df_sequences["ref_pos_start"].astype(int)
    df_sequences["ref_pos_end"] = df_sequences["ref_pos_end"].astype(int)
    df = df.merge(df_sequences, on=["contig_name","ref_pos_start","ref_pos_end","orientation"])

    # fetch rare case: multiple point mutations found in the same sequence results
    # in duplicated rows when merged by pos-on-contig only
    df = df.drop_duplicates(subset=["contig_name","ref_pos_start","ref_pos_end","orientation","sequence"])

    return df


def read_mlst_results(input_file: str) -> pd.DataFrame:
    """
    """
    df = pd.read_csv(input_file, sep=",", header=None, na_values=["-"])
    df["allele_types"] = df.iloc[:,3:].apply(", ".join, axis=1)
    df = df.rename(columns={1: "scheme_name", 2:"sequence_type"})
    return df[["scheme_name","sequence_type","allele_types"]]
