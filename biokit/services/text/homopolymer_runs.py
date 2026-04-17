from itertools import groupby
from typing import Any

from .base import Text
from ...helpers.files import read_and_parse_fasta_seqio


class HomopolymerRuns(Text):
    def __init__(self, args: Any) -> None:
        self.per_base: bool = bool(getattr(args, "per_base", False))
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        if self.fasta is None:
            raise ValueError("fasta cannot be None")
        output_format = self.normalize_output_format(self.output_format)
        records = read_and_parse_fasta_seqio(self.fasta)

        if self.per_base:
            rows = self._per_base_rows(records)
            field_order = ["id", "A", "C", "G", "T"]
            if output_format == "tsv":
                print("\t".join(field_order))
                for row in rows:
                    print("\t".join(str(row[f]) for f in field_order))
                return
            print(self.format_rows(rows, output_format, field_order=field_order))
            return

        rows = self._longest_rows(records)
        field_order = ["id", "length", "base", "position"]
        if output_format == "tsv":
            print("\t".join(field_order))
            for row in rows:
                print("\t".join(str(row[f]) for f in field_order))
            return
        print(self.format_rows(rows, output_format, field_order=field_order))

    @staticmethod
    def _longest_run(seq: str) -> tuple[int, str, int]:
        longest_len = 0
        longest_base = ""
        longest_pos = 0
        pos = 0
        for base, group in groupby(seq):
            length = sum(1 for _ in group)
            if length > longest_len:
                longest_len = length
                longest_base = base
                longest_pos = pos + 1
            pos += length
        return longest_len, longest_base, longest_pos

    @staticmethod
    def _longest_per_base(seq: str) -> dict[str, int]:
        per_base = {"A": 0, "C": 0, "G": 0, "T": 0}
        for base, group in groupby(seq):
            if base in per_base:
                length = sum(1 for _ in group)
                if length > per_base[base]:
                    per_base[base] = length
        return per_base

    def _longest_rows(self, records: Any) -> list[dict[str, Any]]:
        rows = []
        for record in records:
            seq = str(record.seq).upper()
            length, base, position = self._longest_run(seq)
            rows.append(
                {
                    "id": record.id,
                    "length": length,
                    "base": base,
                    "position": position,
                }
            )
        return rows

    def _per_base_rows(self, records: Any) -> list[dict[str, Any]]:
        rows = []
        for record in records:
            seq = str(record.seq).upper()
            per_base = self._longest_per_base(seq)
            row: dict[str, Any] = {"id": record.id}
            row.update(per_base)
            rows.append(row)
        return rows

    def process_args(self, args: Any) -> dict[str, str | None]:
        return dict(fasta=args.fasta, output_format=getattr(args, "format", None))
