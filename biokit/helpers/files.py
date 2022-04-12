import logging
import sys

from Bio import AlignIO, SeqIO

log = logging.getLogger(__name__)


def read_alignment_alignio(fasta):
    try:
        if fasta == '-':
            alignment = AlignIO.read(sys.stdin, "fasta")
        else:
            alignment = AlignIO.read(open(fasta), "fasta")
        return alignment
    except FileNotFoundError:
        raise Exception("Input file could not be read")


def read_and_parse_fasta_seqio(fasta):
    try:
        if fasta == '-':
            records = SeqIO.parse(sys.stdin, "fasta")
        else:
            records = SeqIO.parse(open(fasta), "fasta")
        return records
    except FileNotFoundError:
        raise Exception("Input file could not be read")
