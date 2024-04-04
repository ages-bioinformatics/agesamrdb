#
from Bio import SeqRecord
import zlib

ADAPT_COLS = ["Start","Stop","isBegin","isEnd","start1","end1","start2","end2","orfBegin","orfEnd",
        "start","stop","start_attL","end_attL","start_attR","end_attR","ref_pos_start","ref_pos_end"]

def calc_sequence_hash(sequence: str):
    return hex(zlib.crc32(sequence.encode('ascii')))


def gene_quality_control(seqrecord, ignore_missing_stop=True):
    """
    some basic quality control about the sequence :start and stopcodon, 
    checking wheter we see a frameshift
    """
    comment = []
    if not ignore_missing_stop:
        if not str(seqrecord.seq[-3:]).upper() in ["TAA","TAG","TGA"]:
            comment.append("Missing STOP codon")
    if not str(seqrecord.seq[:3]).upper() in ["ATG", "GTG", "TTG"]:
        comment.append("Missing START codon")
    if not len(seqrecord.seq.replace("-","")) % 3 == 0:
        comment.append("Consensus length not in codon tripplets")

    if len(comment) > 0:
        return ", ".join(comment)
    else:
        return None


def get_or_create(session, model, **kwargs):
    """
    generic function similar to what's known from django ORM
    creates new completely black item if all kwargs are None (assumed that this
    behaviour is allowed in database and as cli-args)
    """
    instance = session.query(model).filter_by(**kwargs).first()
    if not all(v is None for v in kwargs.values()) and instance is not None:
        return instance
    else:
        instance = model(**kwargs)
        session.add(instance)
        session.commit()
        return instance


def apply_offset_to_partial_contigs(df, contig_id_col="seqID"):
    df[[f"new_{contig_id_col}", "offset"]] = df[contig_id_col].str.extract(r"(\S+)\:(\d+)-\d+$", expand=True)
    # check if not applicable (nothing got extracted)
    if not df["offset"].isna().all():
        df["offset"] = df["offset"].fillna(0)
        df["offset"] = df["offset"].astype(int)
        for col in df.columns:
            if col in ADAPT_COLS:
                df[col] = df[col] + df["offset"]
        df[contig_id_col] = df[f"new_{contig_id_col}"]

    df = df.drop(["offset",f"new_{contig_id_col}"], axis=1)
    return df

