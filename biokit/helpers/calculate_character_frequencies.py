from collections import Counter

from ..helpers.files import read_and_parse_fasta_seqio


def calculate_character_frequencies(fasta: str) -> dict:
    """
    determine the frequency of characters in a
    FASTA file
    """
    # get contig lengths
    records = read_and_parse_fasta_seqio(fasta)
    seqs = []
    for record in records:
        if isinstance(record.seq._data.upper(), str):
            seqs.append(record.seq._data.upper())
        else:
            seqs.append(record.seq._data.decode("utf-8").upper())
    res = dict(Counter("".join(seqs)))

    return res
