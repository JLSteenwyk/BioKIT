import sys

from Bio import SeqIO

from .base import Genome


class NumberOfScaffolds(Genome):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        print(self.number_of_scaffolds())

    def process_args(self, args):
        return dict(fasta=args.fasta)
