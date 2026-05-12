from collections import Counter
from typing import Any

from .base import Text
from ...helpers.files import read_and_parse_fasta_seqio


NUCLEOTIDES = ("A", "C", "G", "T")
DINUCLEOTIDES = tuple(a + b for a in NUCLEOTIDES for b in NUCLEOTIDES)


class DinucleotideOdds(Text):
    def __init__(self, args: Any) -> None:
        self.verbose: bool = bool(getattr(args, "verbose", False))
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        if self.fasta is None:
            raise ValueError("fasta cannot be None")
        output_format = self.normalize_output_format(self.output_format)
        records = read_and_parse_fasta_seqio(self.fasta)

        if self.verbose:
            rows = []
            for record in records:
                ratios = self._compute_ratios(str(record.seq))
                row: dict[str, Any] = {"id": record.id}
                for dinuc in DINUCLEOTIDES:
                    row[dinuc] = ratios[dinuc]
                rows.append(row)
            field_order = ["id"] + list(DINUCLEOTIDES)
            if output_format == "tsv":
                print("\t".join(field_order))
                for row in rows:
                    print("\t".join(str(row[f]) for f in field_order))
                return
            print(self.format_rows(rows, output_format, field_order=field_order))
            return

        mono_counts: Counter[str] = Counter()
        di_counts: Counter[str] = Counter()
        for record in records:
            self._accumulate_counts(str(record.seq), mono_counts, di_counts)
        ratios = self._ratios_from_counts(mono_counts, di_counts)

        if output_format == "tsv":
            print("dinucleotide\tO_E")
            for dinuc in DINUCLEOTIDES:
                print(f"{dinuc}\t{ratios[dinuc]}")
            return
        rows = [{"dinucleotide": dn, "O_E": ratios[dn]} for dn in DINUCLEOTIDES]
        print(self.format_rows(rows, output_format, field_order=["dinucleotide", "O_E"]))

    @staticmethod
    def _accumulate_counts(
        sequence: str,
        mono_counts: Counter[str],
        di_counts: Counter[str],
    ) -> None:
        seq = sequence.upper()
        valid = set(NUCLEOTIDES)
        for base in seq:
            if base in valid:
                mono_counts[base] += 1
        for i in range(len(seq) - 1):
            pair = seq[i:i + 2]
            if pair[0] in valid and pair[1] in valid:
                di_counts[pair] += 1

    @classmethod
    def _ratios_from_counts(
        cls,
        mono_counts: Counter[str],
        di_counts: Counter[str],
    ) -> dict[str, float | None]:
        total_mono = sum(mono_counts.values())
        total_di = sum(di_counts.values())
        ratios: dict[str, float | None] = {}
        if total_mono == 0 or total_di == 0:
            for dn in DINUCLEOTIDES:
                ratios[dn] = None
            return ratios
        mono_freq = {b: mono_counts[b] / total_mono for b in NUCLEOTIDES}
        for dn in DINUCLEOTIDES:
            obs = di_counts.get(dn, 0) / total_di
            expected = mono_freq[dn[0]] * mono_freq[dn[1]]
            if expected == 0:
                ratios[dn] = None
            else:
                ratios[dn] = round(obs / expected, 4)
        return ratios

    @classmethod
    def _compute_ratios(cls, sequence: str) -> dict[str, float | None]:
        mono_counts: Counter[str] = Counter()
        di_counts: Counter[str] = Counter()
        cls._accumulate_counts(sequence, mono_counts, di_counts)
        return cls._ratios_from_counts(mono_counts, di_counts)

    def process_args(self, args: Any) -> dict[str, str | None]:
        return dict(fasta=args.fasta, output_format=getattr(args, "format", None))
