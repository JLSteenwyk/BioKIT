from typing import Any

from Bio.SeqUtils.ProtParam import ProteinAnalysis

from .base import Text
from ...helpers.files import read_and_parse_fasta_seqio


class ProteinCharge(Text):
    def __init__(self, args: Any) -> None:
        self.pH: float = float(args.pH) if args.pH is not None else 7.0
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        if self.fasta is None:
            raise ValueError("fasta cannot be None")
        output_format = self.normalize_output_format(self.output_format)
        records = read_and_parse_fasta_seqio(self.fasta)
        if output_format == "tsv":
            for seq_record in records:
                seq_str = str(seq_record.seq).upper()
                analysis = ProteinAnalysis(seq_str)
                charge = round(analysis.charge_at_pH(self.pH), 4)
                print(f"{seq_record.id}\t{charge}")
            return

        rows = []
        for seq_record in records:
            seq_str = str(seq_record.seq).upper()
            analysis = ProteinAnalysis(seq_str)
            charge = round(analysis.charge_at_pH(self.pH), 4)
            rows.append({"id": seq_record.id, "charge": charge, "pH": self.pH})
        rows.sort(key=lambda row: row["id"])
        print(self.format_rows(rows, output_format))

    def process_args(self, args: Any) -> dict[str, str | None]:
        return dict(
            fasta=args.fasta,
            output_format=getattr(args, "format", None),
        )
