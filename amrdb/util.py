#
from Bio import SeqRecord
import zlib

def calc_sequence_hash(sequence: str):
    return hex(zlib.crc32(sequence.encode('ascii')))


def gene_quality_control(seqrecord):
    """
    some basic quality control about the sequence (start, stopcodon, frameshift)
    """
    comment = []
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

