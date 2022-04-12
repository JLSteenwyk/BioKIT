from Bio import SeqIO

from .base import Text
from ...helpers.files import read_and_parse_fasta_seqio


class RemoveFastaEntry(Text):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        record_dict = read_and_parse_fasta_seqio(self.fasta)

        out_records = []

        for i in record_dict:
            if i.name != self.entry:
                out_records.append(i)

        SeqIO.write(out_records, self.output, "fasta")

    def process_args(self, args):
        if args.output is None:
            output = args.fasta + ".pruned.fa"
        else:
            output = args.output

        if output == "-.pruned.fa":
            output = "pruned.fa"

        return dict(fasta=args.fasta, entry=args.entry, output=output)
