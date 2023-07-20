import pandas as pd
from .models import Sample, ISEScanResult, BaktaResult
from .io import read_isescan_results, read_bakta_results, read_resfinder_results, read_pointfinder_results
from .insert import insert_generic_contig_results, insert_into_resfinder_results, insert_into_pointfinder_results, \
        add_new_sequences, add_contig_info

BAKTA_DB_COLUMNS=["product_type","ref_pos_start","ref_pos_end","orientation","gene_name",
        "product","db_xref"]
ISESCAN_DB_COLUMNS = ["family","cluster","is_start_pos","is_end_pos","is_copy_number",
        "ir_start_pos1","ir_end_pos1","ir_start_pos2","ir_end_pos2","score","irId",
        "irLen","nGaps","orf_start_pos","orf_end_pos","orientation","e_value","complete","ov","tir"]


def read_result(inputpath: str, method: str, **kwargs) -> pd.DataFrame:
    """
    interface function that can be universally used to read data into
    a pandas dataframe that contains columns which will be written to database
    params:
    inputpath: string, path to tabular file or output directory of tool
    method: must be string from set: "isescan","bakta","resfinder","pointfinder"
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
    else:
        raise LookupError(f"Method not implemented: {method}")


def insert_into_db(df: pd.DataFrame, method: str, associated_sample: Sample, session: object, **kwargs) -> None:
    """
    interface function to import data of generic pandas.DataFrame format to database
    """
    if method == "isescan":
        return insert_generic_contig_results(df, associated_sample, session, to_db_columns=ISESCAN_DB_COLUMNS, model=ISEScanResult)
    elif method == "bakta":
        return insert_generic_contig_results(df, associated_sample, session, to_db_columns=BAKTA_DB_COLUMNS, model=BaktaResult)
    elif method == "resfinder":
        add_new_sequences(df, session)
        if kwargs.get("assembly_path"):
            df = add_contig_info(df, kwargs["assembly_path"])
        return insert_into_resfinder_results(df, associated_sample, session)
    elif method == "pointfinder":
        return insert_into_pointfinder_results(df, associated_sample, session)
    else:
        raise LookupError (f"Method not implemented: {method}")
