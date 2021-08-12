import sys

from Bio import SeqIO

from .base import Text

class L90(Text):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        # get contig lengths
        records = SeqIO.parse(self.fasta, "fasta")
        contig_lens = []
        for seq_record in records:
            contig_lens.append(len(seq_record))
        
        # sort and reverse contig lengths
        contig_lens = sorted(contig_lens)
        contig_lens.reverse()

        # calculate N50
        sum_contig_lens = sum(contig_lens)
        threshold = sum_contig_lens*.90
        n90 = 0
        l90 = 0

        for contig_len in contig_lens:
            n90 += contig_len
            l90 += 1
            if n90 >= threshold:
                print(l90)
                break
    
    def process_args(self, args):
        return dict(fasta=args.fasta)