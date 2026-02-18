import logging
import sys
from typing import Any, Iterator

log = logging.getLogger(__name__)

_AlignIO: Any | None = None
try:
    from Bio import AlignIO as _AlignIO
except ImportError:  # pragma: no cover - dependency is expected in runtime env
    pass

AlignIO: Any | None = _AlignIO


def read_alignment_alignio(fasta: str) -> Any:
    try:
        if AlignIO is None:
            raise Exception("Biopython AlignIO is not available")

        if fasta == '-':
            alignment = AlignIO.read(sys.stdin, "fasta")
        else:
            with open(fasta) as fasta_handle:
                alignment = AlignIO.read(fasta_handle, "fasta")
        return alignment
    except FileNotFoundError:
        raise Exception("Input file could not be read")


def read_and_parse_fasta_seqio(fasta: str) -> Any:
    try:
        from Bio import SeqIO

        if fasta == '-':
            records = SeqIO.parse(sys.stdin, "fasta")
        else:
            records = SeqIO.parse(fasta, "fasta")
        return records
    except FileNotFoundError:
        raise Exception("Input file could not be read")


def iter_fasta_sequences(fasta: str) -> Iterator[str]:
    """
    iterate through sequence strings in a FASTA file using
    Biopython's lightweight FASTA parser.
    """
    try:
        from Bio.SeqIO.FastaIO import SimpleFastaParser

        if fasta == "-":
            for _, seq in SimpleFastaParser(sys.stdin):
                yield seq
        else:
            with open(fasta) as fasta_handle:
                for _, seq in SimpleFastaParser(fasta_handle):
                    yield seq
    except FileNotFoundError:
        raise Exception("Input file could not be read")


def iter_fasta_entries(fasta: str) -> Iterator[tuple[str, str]]:
    """
    iterate through FASTA entry ids and sequence strings using a
    lightweight parser.
    """
    try:
        from Bio.SeqIO.FastaIO import SimpleFastaParser

        if fasta == "-":
            for title, seq in SimpleFastaParser(sys.stdin):
                yield title.split(None, 1)[0], seq
        else:
            with open(fasta) as fasta_handle:
                for title, seq in SimpleFastaParser(fasta_handle):
                    yield title.split(None, 1)[0], seq
    except FileNotFoundError:
        raise Exception("Input file could not be read")
