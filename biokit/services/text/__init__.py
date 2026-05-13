import importlib
from typing import Any

__all__ = [
    "CharacterFrequency",
    "DinucleotideOdds",
    "Faidx",
    "FastaDeduplication",
    "FileFormatConverter",
    "FindOrfs",
    "GenbankToFasta",
    "HomopolymerRuns",
    "KmerFrequency",
    "MeltingTemperature",
    "MultipleLineToSingleLineFasta",
    "ProteinCharge",
    "ProteinProperties",
    "RenameFastaEntries",
    "RestrictionSites",
    "SampleSequences",
    "ShuffleSequences",
    "RemoveFastaEntry",
    "RemoveShortSequences",
    "ReorderBySequenceLength",
    "SequenceComplement",
    "SequenceLength",
    "SingleLineToMultipleLineFasta",
]

_LAZY_IMPORTS = {
    "CharacterFrequency": ".character_frequency",
    "DinucleotideOdds": ".dinucleotide_odds",
    "Faidx": ".faidx",
    "FastaDeduplication": ".fasta_deduplication",
    "FileFormatConverter": ".file_format_converter",
    "FindOrfs": ".find_orfs",
    "GenbankToFasta": ".genbank_to_fasta",
    "HomopolymerRuns": ".homopolymer_runs",
    "KmerFrequency": ".kmer_frequency",
    "MeltingTemperature": ".melting_temperature",
    "MultipleLineToSingleLineFasta": ".multiple_line_to_single_line_fasta",
    "ProteinCharge": ".protein_charge",
    "ProteinProperties": ".protein_properties",
    "RenameFastaEntries": ".rename_fasta_entries",
    "RestrictionSites": ".restriction_sites",
    "SampleSequences": ".sample_sequences",
    "ShuffleSequences": ".shuffle_sequences",
    "RemoveFastaEntry": ".remove_fasta_entry",
    "RemoveShortSequences": ".remove_short_sequences",
    "ReorderBySequenceLength": ".reorder_by_sequence_length",
    "SequenceComplement": ".sequence_complement",
    "SequenceLength": ".sequence_length",
    "SingleLineToMultipleLineFasta": ".single_line_to_multiple_line_fasta",
}


def __getattr__(name: str) -> Any:
    if name in _LAZY_IMPORTS:
        module = importlib.import_module(_LAZY_IMPORTS[name], __name__)
        return getattr(module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
