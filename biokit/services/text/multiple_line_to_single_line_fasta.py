from Bio import SeqIO

from .base import Text


class MultipleLineToSingleLineFasta(Text):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        # create biopython object of sequences
        records = SeqIO.parse(self.fasta, "fasta")

        res_records = []
        for record in records:
            res_records.append(">" + record.id)
            res_records.append(record.seq._data)

        print('\n'.join(res_records))

    def process_args(self, args):
        return dict(
            fasta=args.fasta
        )
