#
from Bio import SeqRecord
import zlib

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

