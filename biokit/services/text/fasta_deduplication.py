from typing import Any

from Bio import SeqIO

from .base import Text
from ...helpers.files import read_and_parse_fasta_seqio


class FastaDeduplication(Text):
    def __init__(self, args: Any) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        if self.fasta is None or self.output is None:
            raise ValueError("fasta and output cannot be None")
        records = read_and_parse_fasta_seqio(self.fasta)

        seen: set[str] = set()
        unique_records: list[Any] = []

        for record in records:
            seq_upper = str(record.seq).upper()
            if seq_upper not in seen:
                seen.add(seq_upper)
                unique_records.append(record)

        SeqIO.write(unique_records, self.output, "fasta")

    def process_args(self, args: Any) -> dict[str, Any]:
        if args.output is None:
            output = args.fasta + ".dedup.fa"
        else:
            output = args.output

        return dict(
            fasta=args.fasta,
            output=output,
        )
