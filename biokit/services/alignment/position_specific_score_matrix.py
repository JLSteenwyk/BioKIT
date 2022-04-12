from Bio.Align import AlignInfo

from .base import Alignment
from ...helpers.files import read_alignment_alignio


class PositionSpecificScoreMatrix(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment = read_alignment_alignio(self.fasta)
        for seqrecord in alignment:
            seqrecord.seq._data = seqrecord.seq._data.upper()
            seqrecord.seq._data = seqrecord.seq._data.replace("?", "-")
        summary_align = AlignInfo.SummaryInfo(alignment)
        consensus = summary_align.dumb_consensus()
        my_pssm = summary_align.pos_specific_score_matrix(
            consensus, chars_to_ignore=self.ambiguous_character
        )
        print(vars(my_pssm))

    def process_args(self, args):
        return dict(
            fasta=args.fasta,
            ambiguous_character=args.ambiguous_character,
        )
