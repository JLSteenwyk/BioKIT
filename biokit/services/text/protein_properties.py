from typing import Any

from Bio.SeqUtils.ProtParam import ProteinAnalysis

from .base import Text
from ...helpers.files import read_and_parse_fasta_seqio


PROPERTY_FIELDS = (
    "id",
    "length",
    "molecular_weight",
    "isoelectric_point",
    "gravy",
    "aromaticity",
    "instability_index",
)


class ProteinProperties(Text):
    def __init__(self, args: Any) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        if self.fasta is None:
            raise ValueError("fasta cannot be None")
        output_format = self.normalize_output_format(self.output_format)
        records = read_and_parse_fasta_seqio(self.fasta)

        rows = []
        for seq_record in records:
            seq_str = self._clean_sequence(str(seq_record.seq))
            row: dict[str, Any] = {"id": seq_record.id, "length": len(seq_str)}
            row.update(self._compute_properties(seq_str))
            rows.append(row)

        if output_format == "tsv":
            print("\t".join(PROPERTY_FIELDS))
            for row in rows:
                print("\t".join(str(row[f]) for f in PROPERTY_FIELDS))
            return
        print(self.format_rows(rows, output_format, field_order=list(PROPERTY_FIELDS)))

    @staticmethod
    def _clean_sequence(seq: str) -> str:
        upper = seq.upper()
        return "".join(c for c in upper if c not in {"*", "-"})

    @staticmethod
    def _compute_properties(seq_str: str) -> dict[str, Any]:
        if not seq_str:
            return {
                "molecular_weight": None,
                "isoelectric_point": None,
                "gravy": None,
                "aromaticity": None,
                "instability_index": None,
            }
        analysis = ProteinAnalysis(seq_str)
        try:
            mw = round(analysis.molecular_weight(), 4)
        except (KeyError, ValueError):
            mw = None
        try:
            pi = round(analysis.isoelectric_point(), 4)
        except (KeyError, ValueError):
            pi = None
        try:
            gravy = round(analysis.gravy(), 4)
        except (KeyError, ValueError):
            gravy = None
        try:
            aromaticity = round(analysis.aromaticity(), 4)
        except (KeyError, ValueError):
            aromaticity = None
        try:
            instability = round(analysis.instability_index(), 4)
        except (KeyError, ValueError):
            instability = None
        return {
            "molecular_weight": mw,
            "isoelectric_point": pi,
            "gravy": gravy,
            "aromaticity": aromaticity,
            "instability_index": instability,
        }

    def process_args(self, args: Any) -> dict[str, str | None]:
        return dict(
            fasta=args.fasta,
            output_format=getattr(args, "format", None),
        )
