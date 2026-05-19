from itertools import combinations
from typing import Any

from .base import Alignment
from ...helpers.files import read_alignment_alignio


VALID_BASES = frozenset("ACGT")
TRANSITION_PAIRS = frozenset({
    ("A", "G"), ("G", "A"),
    ("C", "T"), ("T", "C"),
})


class TransitionTransversionRatio(Alignment):
    def __init__(self, args: Any) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        if self.fasta is None:
            raise ValueError("fasta cannot be None")
        output_format = self.normalize_output_format(self.output_format)
        alignment = read_alignment_alignio(self.fasta)
        aln_len = alignment.get_alignment_length()

        transitions, transversions, site_classifications = self._scan_alignment(
            alignment, aln_len
        )
        ratio = self._compute_ratio(transitions, transversions)

        if self.verbose:
            self._emit_verbose(site_classifications, output_format)
            return

        self._emit_summary(transitions, transversions, ratio, output_format)

    @classmethod
    def _scan_alignment(
        cls, alignment: Any, aln_len: int
    ) -> tuple[int, int, list[tuple[int, str]]]:
        transitions = 0
        transversions = 0
        site_classifications: list[tuple[int, str]] = []

        for i in range(aln_len):
            column = alignment[:, i].upper()
            site_ti, site_tv, label = cls._classify_column(column)
            transitions += site_ti
            transversions += site_tv
            site_classifications.append((i + 1, label))

        return transitions, transversions, site_classifications

    @classmethod
    def _classify_column(cls, column: str) -> tuple[int, int, str]:
        if any(base in ("-", "?") for base in column):
            return 0, 0, "gap"

        non_gap_bases = [base for base in column if base in VALID_BASES]
        unique_bases = set(non_gap_bases)
        if len(unique_bases) <= 1:
            return 0, 0, "constant"

        site_ti = 0
        site_tv = 0
        has_tv = False
        for a, b in combinations(non_gap_bases, 2):
            if a == b:
                continue
            if (a, b) in TRANSITION_PAIRS:
                site_ti += 1
            else:
                site_tv += 1
                has_tv = True

        label = "transversion" if has_tv else "transition"
        return site_ti, site_tv, label

    @staticmethod
    def _compute_ratio(transitions: int, transversions: int) -> float | None:
        if transversions == 0:
            return None
        return round(transitions / transversions, 4)

    def _emit_summary(
        self, transitions: int, transversions: int,
        ratio: float | None, output_format: str,
    ) -> None:
        if output_format == "tsv":
            print(f"{transitions}\t{transversions}\t{ratio}")
            return
        print(self.format_object(
            {
                "transitions": transitions,
                "transversions": transversions,
                "ratio": ratio,
            },
            output_format,
        ))

    def _emit_verbose(
        self,
        site_classifications: list[tuple[int, str]],
        output_format: str,
    ) -> None:
        if output_format == "tsv":
            for position, label in site_classifications:
                print(f"{position}\t{label}")
            return
        rows = [
            {"position": position, "classification": label}
            for position, label in site_classifications
        ]
        print(self.format_rows(
            rows, output_format, field_order=["position", "classification"]
        ))

    def process_args(self, args: Any) -> dict[str, Any | None]:
        return dict(
            fasta=args.fasta,
            verbose=args.verbose,
            output_format=getattr(args, "format", None),
        )
