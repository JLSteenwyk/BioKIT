from Bio import SeqIO

from .base import Text
from ...helpers.files import read_and_parse_fasta_seqio


class RemoveShortSequences(Text):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        record_dict = read_and_parse_fasta_seqio(self.fasta)

        out_records = []

        for i in record_dict:
            if len(i.seq) > self.threshold:
                out_records.append(i)

        SeqIO.write(out_records, self.output, "fasta")

    def process_args(self, args):
        if args.threshold is None:
            threshold = 500
        else:
            threshold = int(args.threshold)

        if args.output is None:
            output = args.fasta + ".long_seqs.fa"
        else:
            output = args.output

        return dict(
            fasta=args.fasta,
            output=output,
            threshold=threshold
        )
