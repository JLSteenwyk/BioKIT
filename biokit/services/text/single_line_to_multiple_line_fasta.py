from .base import Text
from ...helpers.files import read_and_parse_fasta_seqio


class SingleLineToMultipleLineFasta(Text):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        # create biopython object of sequences
        records = read_and_parse_fasta_seqio(self.fasta)

        res_records = []
        for record in records:
            res_records.append(">" + record.id)
            for i in range(0, len(record.seq._data), 60):
                res_records.append(record.seq._data[i:i + 60])

        print('\n'.join(res_records))

    def process_args(self, args):
        return dict(
            fasta=args.fasta
        )
