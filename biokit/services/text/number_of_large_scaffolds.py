import sys

from Bio import SeqIO

from .base import Text

class NumberOfLargeScaffolds(Text):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        records = SeqIO.parse(self.fasta, "fasta")
        cnt = 0
        total_len = 0
        for seq_record in records:
            seq_len = len(seq_record)
            if seq_len > self.threshold:
                cnt += 1
                total_len += seq_len

        print(f"{cnt}\t{total_len}")
    
    def process_args(self, args):
        if args.threshold is None:
            threshold = 500
        else:
            threshold = int(args.threshold)

        return dict(
            fasta=args.fasta,
            threshold=threshold,
        )