import importlib

__all__ = [
    "GCContentFirstPosition",
    "GCContentSecondPosition",
    "GCContentThirdPosition",
    "GeneWiseRelativeSynonymousCodonUsage",
    "RelativeSynonymousCodonUsage",
    "TranslateSequence",
]

_LAZY_IMPORTS = {
    "GCContentFirstPosition": ".gc_content_first_position",
    "GCContentSecondPosition": ".gc_content_second_position",
    "GCContentThirdPosition": ".gc_content_third_position",
    "GeneWiseRelativeSynonymousCodonUsage": ".gene_wise_relative_synonymous_codon_usage",
    "RelativeSynonymousCodonUsage": ".relative_synonymous_codon_usage",
    "TranslateSequence": ".translate_sequence",
}


def __getattr__(name):
    if name in _LAZY_IMPORTS:
        module = importlib.import_module(_LAZY_IMPORTS[name], __name__)
        return getattr(module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
