import importlib

__all__ = [
    "FastQReadLengths",
    "SubsetPEFastQReads",
    "SubsetSEFastQReads",
    "TrimPEAdaptersFastQ",
    "TrimPEFastQ",
    "TrimSEAdaptersFastQ",
    "TrimSEFastQ",
]

_LAZY_IMPORTS = {
    "FastQReadLengths": ".fastq_read_lengths",
    "SubsetPEFastQReads": ".subset_pe_fastq_reads",
    "SubsetSEFastQReads": ".subset_se_fastq_reads",
    "TrimPEAdaptersFastQ": ".trim_pe_adapters_fastq",
    "TrimPEFastQ": ".trim_pe_fastq",
    "TrimSEAdaptersFastQ": ".trim_se_adapters_fastq",
    "TrimSEFastQ": ".trim_se_fastq",
}


def __getattr__(name):
    if name in _LAZY_IMPORTS:
        module = importlib.import_module(_LAZY_IMPORTS[name], __name__)
        return getattr(module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
