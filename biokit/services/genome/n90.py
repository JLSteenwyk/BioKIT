import sys

from Bio import SeqIO

from .base import Genome


class N90(Genome):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        print(self.calc_n90())

    def process_args(self, args):
        return dict(fasta=args.fasta)
