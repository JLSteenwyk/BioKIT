import sys

from Bio import SeqIO

from .base import Genome


class LongestScaffold(Genome):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        # get contig lengths
        records = SeqIO.parse(self.fasta, "fasta")
        contig_lens = []
        for seq_record in records:
            contig_lens.append(len(seq_record))

        # get longest contig length
        print(max(contig_lens))

    def process_args(self, args):
        return dict(fasta=args.fasta)
