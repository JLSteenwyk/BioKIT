import re

from Bio import AlignIO
from Bio.Align import AlignInfo

from .base import Alignment


class ConsensusSequence(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment = AlignIO.read(self.fasta, "fasta")
        summary_align = AlignInfo.SummaryInfo(alignment)

        if not self.ambiguous_character:
            ambiguous_character = "N"
        else:
            ambiguous_character = self.ambiguous_character

        if not self.threshold:
            threshold = 0.7
        else:
            threshold = float(self.threshold)

        consensus = summary_align.dumb_consensus(
            threshold=threshold, ambiguous=ambiguous_character
        )

        header = ">" + re.sub("^.*/", "", str(self.fasta)) + ".consensus"
        print(f"{header}\n{consensus}")

    def process_args(self, args):
        return dict(
            fasta=args.fasta,
            threshold=args.threshold,
            ambiguous_character=args.ambiguous_character,
        )
