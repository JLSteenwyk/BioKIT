import importlib

__all__ = [
    "CharacterFrequency",
    "Faidx",
    "FileFormatConverter",
    "MultipleLineToSingleLineFasta",
    "RenameFastaEntries",
    "RemoveFastaEntry",
    "RemoveShortSequences",
    "ReorderBySequenceLength",
    "SequenceComplement",
    "SequenceLength",
    "SingleLineToMultipleLineFasta",
]

_LAZY_IMPORTS = {
    "CharacterFrequency": ".character_frequency",
    "Faidx": ".faidx",
    "FileFormatConverter": ".file_format_converter",
    "MultipleLineToSingleLineFasta": ".multiple_line_to_single_line_fasta",
    "RenameFastaEntries": ".rename_fasta_entries",
    "RemoveFastaEntry": ".remove_fasta_entry",
    "RemoveShortSequences": ".remove_short_sequences",
    "ReorderBySequenceLength": ".reorder_by_sequence_length",
    "SequenceComplement": ".sequence_complement",
    "SequenceLength": ".sequence_length",
    "SingleLineToMultipleLineFasta": ".single_line_to_multiple_line_fasta",
}


def __getattr__(name):
    if name in _LAZY_IMPORTS:
        module = importlib.import_module(_LAZY_IMPORTS[name], __name__)
        return getattr(module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
