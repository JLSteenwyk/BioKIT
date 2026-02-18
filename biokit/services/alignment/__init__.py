import importlib

__all__ = [
    "AlignmentLength",
    "AlignmentRecoding",
    "AlignmentSummary",
    "ConsensusSequence",
    "ConstantSites",
    "ParsimonyInformativeSites",
    "PositionSpecificScoreMatrix",
    "VariableSites",
]

_LAZY_IMPORTS = {
    "AlignmentLength": ".alignment_length",
    "AlignmentRecoding": ".alignment_recoding",
    "AlignmentSummary": ".alignment_summary",
    "ConsensusSequence": ".consensus_sequence",
    "ConstantSites": ".constant_sites",
    "ParsimonyInformativeSites": ".parsimony_informative_sites",
    "PositionSpecificScoreMatrix": ".position_specific_score_matrix",
    "VariableSites": ".variable_sites",
}


def __getattr__(name):
    if name in _LAZY_IMPORTS:
        module = importlib.import_module(_LAZY_IMPORTS[name], __name__)
        return getattr(module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
