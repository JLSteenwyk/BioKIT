from .base import Text

from ...helpers.files import read_and_parse_fasta_seqio


class MultipleLineToSingleLineFasta(Text):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        # create biopython object of sequences
        records = read_and_parse_fasta_seqio(self.fasta)

        res_records = []
        for record in records:
            res_records.append(">" + record.id)
            res_records.append(record.seq._data)

        print('\n'.join(res_records))

    def process_args(self, args):
        return dict(
            fasta=args.fasta
        )
