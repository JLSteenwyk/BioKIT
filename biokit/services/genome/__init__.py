import importlib

__all__ = [
    "GCContent",
    "GenomeAssemblyMetrics",
    "L50",
    "L90",
    "LongestScaffold",
    "N50",
    "N90",
    "NumberOfLargeScaffolds",
    "NumberOfScaffolds",
    "SumOfScaffoldLengths",
]

_LAZY_IMPORTS = {
    "GCContent": ".gc_content",
    "GenomeAssemblyMetrics": ".genome_assembly_metrics",
    "L50": ".l50",
    "L90": ".l90",
    "LongestScaffold": ".longest_scaffold",
    "N50": ".n50",
    "N90": ".n90",
    "NumberOfLargeScaffolds": ".number_of_large_scaffolds",
    "NumberOfScaffolds": ".number_of_scaffolds",
    "SumOfScaffoldLengths": ".sum_of_scaffold_lengths",
}


def __getattr__(name):
    if name in _LAZY_IMPORTS:
        module = importlib.import_module(_LAZY_IMPORTS[name], __name__)
        return getattr(module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
