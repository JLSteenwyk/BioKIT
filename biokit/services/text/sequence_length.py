from typing import Any

from .base import Text
from ...helpers.files import read_and_parse_fasta_seqio


class SequenceLength(Text):
    def __init__(self, args: Any) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        if self.fasta is None:
            raise ValueError("fasta cannot be None")
        output_format = self.normalize_output_format(self.output_format)
        records = read_and_parse_fasta_seqio(self.fasta)
        if output_format == "tsv":
            for seq_record in records:
                print(f"{seq_record.id}\t{len(seq_record)}")
            return

        rows = [
            {"id": seq_record.id, "length": len(seq_record)}
            for seq_record in records
        ]
        rows.sort(key=lambda row: row["id"])
        print(self.format_rows(rows, output_format))

    def process_args(self, args: Any) -> dict[str, str | None]:
        return dict(fasta=args.fasta, output_format=getattr(args, "format", None))
