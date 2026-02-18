from collections import Counter
from typing import Dict

from ..helpers.files import iter_fasta_sequences


def calculate_character_frequencies(fasta: str) -> Dict[str, int]:
    """
    determine the frequency of characters in a
    FASTA file
    """
    counts: Counter[str] = Counter()
    for seq in iter_fasta_sequences(fasta):
        counts.update(seq.upper())
    return dict(counts)
