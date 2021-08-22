import logging

from Bio import AlignIO

log = logging.getLogger(__name__)


def read_alignment_alignio(fasta):
    try:
        alignment = AlignIO.read(open(fasta), "fasta")
        return alignment
    except FileNotFoundError:
        raise Exception("Input file could not be read")
