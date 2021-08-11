import sys

from Bio import SeqIO

from .base import Text

class NumberOfScaffolds(Text):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        # get contig lengths
        records = SeqIO.parse(self.fasta, "fasta")
        cnt = 0
        for record in records:
            cnt += 1
        print(cnt)
    
    def process_args(self, args):
        return dict(fasta=args.fasta)