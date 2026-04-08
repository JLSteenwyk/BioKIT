from typing import Any

from Bio.SeqUtils.MeltingTemp import Tm_NN

from .base import Text
from ...helpers.files import read_and_parse_fasta_seqio


class MeltingTemperature(Text):
    def __init__(self, args: Any) -> None:
        self.na_conc: float = (
            float(args.na) if getattr(args, "na", None) is not None else 50.0
        )
        self.oligo_conc: float = (
            float(args.oligo_conc)
            if getattr(args, "oligo_conc", None) is not None
            else 250.0
        )
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        if self.fasta is None:
            raise ValueError("fasta cannot be None")
        output_format = self.normalize_output_format(self.output_format)
        records = read_and_parse_fasta_seqio(self.fasta)

        if output_format == "tsv":
            for record in records:
                tm = round(
                    Tm_NN(record.seq, Na=self.na_conc, dnac1=self.oligo_conc, dnac2=self.oligo_conc),
                    4,
                )
                print(f"{record.id}\t{tm}")
            return

        rows = []
        for record in records:
            tm = round(
                Tm_NN(record.seq, Na=self.na_conc, dnac1=self.oligo_conc, dnac2=self.oligo_conc),
                4,
            )
            rows.append({"id": record.id, "Tm": tm})
        rows.sort(key=lambda row: row["id"])
        print(self.format_rows(rows, output_format))

    def process_args(self, args: Any) -> dict[str, str | None]:
        return dict(
            fasta=args.fasta,
            output_format=getattr(args, "format", None),
        )
