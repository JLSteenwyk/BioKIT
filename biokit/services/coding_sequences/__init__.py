import importlib
from typing import Any

__all__ = [
    "EffectiveNumberOfCodons",
    "GCContentFirstPosition",
    "GCContentFourFoldDegenerateSites",
    "GCContentSecondPosition",
    "GCContentThirdPosition",
    "GeneWiseRelativeSynonymousCodonUsage",
    "NeutralityPlot",
    "RelativeSynonymousCodonUsage",
    "TranslateSequence",
]

_LAZY_IMPORTS = {
    "EffectiveNumberOfCodons": ".effective_number_of_codons",
    "GCContentFirstPosition": ".gc_content_first_position",
    "GCContentFourFoldDegenerateSites": ".gc_content_four_fold_degenerate_sites",
    "GCContentSecondPosition": ".gc_content_second_position",
    "GCContentThirdPosition": ".gc_content_third_position",
    "GeneWiseRelativeSynonymousCodonUsage": ".gene_wise_relative_synonymous_codon_usage",
    "NeutralityPlot": ".neutrality_plot",
    "RelativeSynonymousCodonUsage": ".relative_synonymous_codon_usage",
    "TranslateSequence": ".translate_sequence",
}


def __getattr__(name: str) -> Any:
    if name in _LAZY_IMPORTS:
        module = importlib.import_module(_LAZY_IMPORTS[name], __name__)
        return getattr(module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
