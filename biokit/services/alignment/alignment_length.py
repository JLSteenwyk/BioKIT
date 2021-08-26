from .base import Alignment
from ...helpers.files import read_alignment_alignio


class AlignmentLength(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment = read_alignment_alignio(self.fasta)

        # alignment length
        print(alignment.get_alignment_length())

    def process_args(self, args):
        return dict(fasta=args.fasta)
