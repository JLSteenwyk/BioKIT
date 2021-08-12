import sys

from Bio import SeqIO

from .base import Text

class N50(Text):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        # get contig lengths
        records = SeqIO.parse(self.fasta, "fasta")
        contig_lens = []
        for seq_record in records:
            contig_lens.append(len(seq_record))
        
        # sort and reverse contig lengths
        contig_lens.sort(reverse=True)

        # calculate N50
        sum_contig_lens = sum(contig_lens)
        n50_threshold = sum_contig_lens*.50
        curr = 0

        for contig_len in contig_lens:
            curr += contig_len
            if curr >= n50_threshold:
                print(contig_len)
                break
    
    def process_args(self, args):
        return dict(fasta=args.fasta)