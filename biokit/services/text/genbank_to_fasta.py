import sys
from typing import Any

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .base import Text


class GenbankToFasta(Text):
    def __init__(self, args: Any) -> None:
        self.feature_types: list[str] | None = self._parse_feature_types(
            getattr(args, "feature_type", None)
        )
        self.translate: bool = bool(getattr(args, "translate", False))
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        if self.input_file is None:
            raise ValueError("input file cannot be None")

        records = list(SeqIO.parse(self.input_file, "genbank"))
        out_records: list[SeqRecord] = []

        if self.feature_types is None:
            for record in records:
                out_records.append(
                    SeqRecord(
                        record.seq,
                        id=record.id,
                        description=record.description or "",
                    )
                )
        else:
            wanted = {ft.lower() for ft in self.feature_types}
            for record in records:
                counters: dict[str, int] = {}
                for feature in record.features:
                    if feature.type.lower() not in wanted:
                        continue
                    counters[feature.type] = counters.get(feature.type, 0) + 1
                    out_records.append(
                        self._build_feature_record(
                            record, feature, counters[feature.type]
                        )
                    )

        self._write_fasta(out_records)

    def _build_feature_record(
        self, record: SeqRecord, feature: Any, index: int
    ) -> SeqRecord:
        quals = feature.qualifiers
        identifier = (
            self._first_qual(quals, "locus_tag")
            or self._first_qual(quals, "gene")
            or self._first_qual(quals, "protein_id")
            or f"{record.id}_{feature.type}_{index}"
        )

        if self.translate and feature.type == "CDS":
            protein = self._first_qual(quals, "translation")
            if protein is not None:
                seq = Seq(protein)
            else:
                seq = feature.extract(record.seq).translate(to_stop=True)
        else:
            seq = feature.extract(record.seq)

        start = int(feature.location.start) + 1
        end = int(feature.location.end)
        strand = "+" if feature.location.strand in (1, None) else "-"
        description = f"{feature.type} {record.id}:{start}-{end}({strand})"

        return SeqRecord(seq, id=identifier, description=description)

    def _write_fasta(self, records: list[SeqRecord]) -> None:
        if self.output_file:
            with open(self.output_file, "w") as handle:
                for record in records:
                    handle.write(self._format_record(record))
            return
        for record in records:
            sys.stdout.write(self._format_record(record))

    @staticmethod
    def _format_record(record: SeqRecord) -> str:
        header = record.id
        if record.description and record.description != record.id:
            header = f"{record.id} {record.description}"
        return f">{header}\n{str(record.seq)}\n"

    @staticmethod
    def _first_qual(qualifiers: dict[str, list[str]], key: str) -> str | None:
        values = qualifiers.get(key)
        if values:
            return values[0]
        return None

    @staticmethod
    def _parse_feature_types(value: str | None) -> list[str] | None:
        if value is None:
            return None
        parts = [part.strip() for part in value.split(",") if part.strip()]
        return parts or None

    def process_args(self, args: Any) -> dict[str, Any]:
        return dict(
            input_file=args.genbank,
            output_file=getattr(args, "output", None),
        )
