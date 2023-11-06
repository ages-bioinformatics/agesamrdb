import pandas as pd
from .models import Sample, ISEScanResult, BaktaResult, MobTyperResult, PlasmidfinderResult
from .io import read_isescan_results, read_bakta_results, read_resfinder_results, \
        read_pointfinder_results, read_mobtyper_results, read_plasmidfinder_results
from .insert import insert_generic_contig_results, insert_into_resfinder_results, \
        insert_into_pointfinder_results, add_new_sequences, add_contig_info


# column names that are actually imported in database as constants
# used for generic "insert_generic_contig_results" function [currently bakta,
# isescan, mobtyper]
BAKTA_DB_COLUMNS=["product_type","ref_pos_start","ref_pos_end","orientation",
        "gene_name","product","db_xref"]
ISESCAN_DB_COLUMNS = ["family","cluster","is_start_pos","is_end_pos",
        "is_copy_number","ir_start_pos1","ir_end_pos1","ir_start_pos2",
        "ir_end_pos2","score","irId","irLen","nGaps","orf_start_pos",
        "orf_end_pos","orientation","e_value","complete","ov","tir"]
MOBTYPER_DB_COLUMNS = ["gc_content","rep_type","rep_type_accession",
        "relaxase_type","relaxase_type_accession","mpf_type",
        "mpf_type_accession","orit_type","orit_accession","predicted_mobility",
        "mash_nearest_neighbor","mash_neighbor_distance",
        "mash_neighbor_identification","primary_cluster_id",
        "secondary_cluster_id","predicted_host_range_overall_rank",
        "predicted_host_range_overall_name","observed_host_range_ncbi_rank",
        "observed_host_range_ncbi_name","reported_host_range_lit_rank",
        "reported_host_range_lit_name","associated_pmid"]
PLASMIDFINDER_DB_COLUMNS = ['database_name', 'plasmid', 'identity', 'note',
        'accession_number', 'query_length', 'template_length', 'ref_pos_start',
        'ref_pos_end']

def read_result(inputpath: str, method: str, **kwargs) -> pd.DataFrame:
    """
    interface function that can be universally used to read data into
    a pandas dataframe that contains columns which will be written to database
    params:
    inputpath: string, path to tabular file or output directory of tool
    method: must be string from set: isescan,bakta,resfinder,pointfinder,
        mobtyper
    kwargs: forwarded to read_resfinder_results
    """
    if method == "isescan":
        return read_isescan_results(inputpath)
    elif method == "bakta":
        return read_bakta_results(inputpath)
    elif method == "resfinder":
        return read_resfinder_results(inputpath, **kwargs)
    elif method == "pointfinder":
        return read_pointfinder_results(inputpath)
    elif method == "mobtyper":
        return read_mobtyper_results(inputpath)
    elif method == "plasmidfinder":
        return read_plasmidfinder_results(inputpath)
    else:
        raise LookupError(f"Method not implemented: {method}")


def insert_into_db(df: pd.DataFrame, method: str, associated_sample: Sample, 
        session: object, **kwargs) -> None:
    """
    interface function to import data of generic pandas.DataFrame format to db
    params:
    df: pandas DataFrame that contains all required columns
    method: must be string from set "isescan","bakta","resfinder","pointfinder"
    associated_sample: sqlalchemy instance of class Sample
    session: sqlalchemy session object
    assembly_path: for resfinder to parse assembly and create contig entries
    """
    if method == "isescan":
        return insert_generic_contig_results(df, associated_sample, session, 
                to_db_columns=ISESCAN_DB_COLUMNS, model=ISEScanResult)
    elif method == "bakta":
        return insert_generic_contig_results(df, associated_sample, session,
                to_db_columns=BAKTA_DB_COLUMNS, model=BaktaResult)
    elif method == "resfinder":
        add_new_sequences(df, session)
        if kwargs.get("assembly_path"):
            df = add_contig_info(df, kwargs["assembly_path"])
        return insert_into_resfinder_results(df, associated_sample, session)
    elif method == "pointfinder":
        return insert_into_pointfinder_results(df, associated_sample, session)
    elif method == "mobtyper":
        return insert_generic_contig_results(df, associated_sample, session,
                to_db_columns=MOBTYPER_DB_COLUMNS, model=MobTyperResult)
    elif method == "plasmidfinder":
    
        return insert_generic_contig_results(df, associated_sample, session,
                to_db_columns=PLASMIDFINDER_DB_COLUMNS, model=PlasmidfinderResult, 
                create_contig=True, contig_name_col='contig_name')
    else:
        raise LookupError (f"Method not implemented: {method}")
