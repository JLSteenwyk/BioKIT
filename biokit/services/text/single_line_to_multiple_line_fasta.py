from typing import Any

from .base import Text
from ...helpers.files import read_and_parse_fasta_seqio


class SingleLineToMultipleLineFasta(Text):
    def __init__(self, args: Any) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        # create biopython object of sequences
        if self.fasta is None:
            raise ValueError("fasta cannot be None")
        records = read_and_parse_fasta_seqio(self.fasta)

        res_records: list[str] = []
        for record in records:
            res_records.append(">" + record.id)
            sequence = str(record.seq)
            for i in range(0, len(sequence), 60):
                res_records.append(sequence[i:i + 60])

        print('\n'.join(res_records))

    def process_args(self, args: Any) -> dict[str, str]:
        return dict(
            fasta=args.fasta
        )
