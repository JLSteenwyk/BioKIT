import sys

from .base import Genome


class N50(Genome):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        print(self.calc_n50())

    def process_args(self, args):
        return dict(fasta=args.fasta)
