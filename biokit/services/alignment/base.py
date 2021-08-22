from ..base import BaseService


class Alignment(BaseService):
    def __init__(
        self,
        *args,
        fasta=None,
        ambiguous_character=None,
        threshold=None,
    ):
        self.fasta = fasta
        self.ambiguous_character = ambiguous_character
        self.threshold = threshold
