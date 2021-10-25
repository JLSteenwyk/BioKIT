from Bio import SeqIO
from collections import Counter


def calculate_character_frequencies(fasta: str) -> dict:
    """
    determine the frequency of characters in a
    FASTA file
    """
    # get contig lengths
    records = SeqIO.parse(fasta, "fasta")
    seqs = []
    for record in records:
        if isinstance(record.seq._data.upper(), str):
            seqs.append(record.seq._data.upper())
        else:
            seqs.append(record.seq._data.decode("utf-8").upper())
    res = dict(Counter("".join(seqs)))

    return res
