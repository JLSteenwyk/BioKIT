import sys

from Bio import SeqIO

from .base import Genome


class SumOfScaffoldLengths(Genome):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        print(self.sum_of_scaffold_lengths())

    def process_args(self, args):
        return dict(fasta=args.fasta)
