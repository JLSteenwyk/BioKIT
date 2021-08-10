import sys

from Bio import SeqIO

from .base import Text

class L50(Text):
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
        threshold = sum_contig_lens*.50
        n50 = 0
        l50 = 0

        for contig_len in contig_lens:
            n50 += contig_len
            l50 += 1
            if n50 >= threshold:
                print(l50)
                break
    
    def process_args(self, args):
        return dict(fasta=args.fasta)