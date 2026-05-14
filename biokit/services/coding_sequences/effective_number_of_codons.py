from collections import Counter, defaultdict
from typing import Any

from .base import CodingSequence
from ...helpers.files import iter_fasta_entries


ENC_MAX = 61.0


class EffectiveNumberOfCodons(CodingSequence):
    def __init__(self, args: Any) -> None:
        self.plot: bool = bool(getattr(args, "plot", False))
        self.output_file: str | None = getattr(args, "output", None)
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        if self.fasta is None:
            raise ValueError("fasta cannot be None")
        output_format = self.normalize_output_format(self.output_format)
        translation_table = self.read_translation_table(self.translation_table)
        synonymous_families = self._synonymous_families(translation_table)

        per_gene = self._compute_per_gene(translation_table, synonymous_families)

        if self.plot:
            self._make_plot(per_gene, synonymous_families)

        if output_format == "tsv":
            for row in per_gene:
                if self.verbose:
                    print(f"{row['gene_id']}\t{row['ENC']}\t{row['GC3']}")
                else:
                    print(f"{row['gene_id']}\t{row['ENC']}")
            return

        if self.verbose:
            field_order = ["gene_id", "ENC", "GC3"]
        else:
            field_order = ["gene_id", "ENC"]
        rows = [{k: row[k] for k in field_order} for row in per_gene]
        print(self.format_rows(rows, output_format, field_order=field_order))

    @staticmethod
    def _synonymous_families(
        translation_table: dict[str, str],
    ) -> dict[str, list[str]]:
        """Map amino acid → list of synonymous codons (excludes stops)."""
        families: dict[str, list[str]] = defaultdict(list)
        for codon, aa in translation_table.items():
            if aa == "*":
                continue
            families[aa].append(codon)
        return dict(families)

    def _compute_per_gene(
        self,
        translation_table: dict[str, str],
        synonymous_families: dict[str, list[str]],
    ) -> list[dict[str, Any]]:
        assert self.fasta is not None
        valid_codons = set(translation_table.keys())
        results: list[dict[str, Any]] = []
        for seq_id, seq in iter_fasta_entries(self.fasta):
            sequence = seq.upper().replace("T", "U")
            if len(sequence) % 3 != 0:
                continue
            codon_counts: Counter[str] = Counter()
            pos3_chars: list[str] = []
            for i in range(0, len(sequence), 3):
                codon = sequence[i:i + 3]
                if codon in valid_codons and translation_table[codon] != "*":
                    codon_counts[codon] += 1
                    pos3_chars.append(codon[2])
            if not codon_counts:
                continue
            enc = self._compute_enc(codon_counts, synonymous_families)
            # GC3 here is computed only over codons that pass our filter
            # (skip non-translatable and stop codons) to stay consistent
            # with the gene's contribution to ENC.
            pos3 = "".join(pos3_chars)
            gc3 = self._gc_fraction(pos3)
            results.append({
                "gene_id": seq_id,
                "ENC": round(enc, 4) if enc is not None else None,
                "GC3": round(gc3, 4),
            })
        return results

    @staticmethod
    def _gc_fraction(seq: str) -> float:
        if not seq:
            return 0.0
        # treat both DNA and RNA G/C
        gc = sum(1 for c in seq if c in ("G", "C"))
        return gc / len(seq)

    @classmethod
    def _compute_enc(
        cls,
        codon_counts: Counter[str],
        synonymous_families: dict[str, list[str]],
    ) -> float | None:
        # Compute F per amino acid (only for families with codons observed)
        f_by_family_size: dict[int, list[float]] = defaultdict(list)
        family_size_counts: dict[int, int] = defaultdict(int)
        for aa, codons in synonymous_families.items():
            family_size = len(codons)
            family_size_counts[family_size] += 1
            if family_size == 1:
                continue
            n = sum(codon_counts.get(c, 0) for c in codons)
            if n == 0:
                continue
            if n == 1:
                f = 1.0
            else:
                squared_sum = sum(
                    (codon_counts.get(c, 0) / n) ** 2 for c in codons
                )
                f = (n * squared_sum - 1) / (n - 1)
            # Clamp F into [1/family_size, 1] (numerical safety)
            f = max(min(f, 1.0), 1.0 / family_size)
            f_by_family_size[family_size].append(f)

        if not f_by_family_size and family_size_counts[1] == 0:
            return None

        mean_f: dict[int, float] = {
            size: sum(values) / len(values) for size, values in f_by_family_size.items()
        }

        # Wright's interpolation: if F̄₃ is missing, estimate from F̄₂ and F̄₄
        # (mean). More generally, for any missing family size present in the
        # genetic code, fall back to a sensible estimate.
        all_sizes = [size for size, count in family_size_counts.items() if count > 0]
        for size in all_sizes:
            if size == 1 or size in mean_f:
                continue
            neighbors = [mean_f[s] for s in mean_f.keys()]
            if not neighbors:
                continue
            mean_f[size] = sum(neighbors) / len(neighbors)

        enc = float(family_size_counts.get(1, 0))
        for size, count in family_size_counts.items():
            if size == 1 or count == 0:
                continue
            f_bar = mean_f.get(size)
            if f_bar is None or f_bar <= 0:
                continue
            enc += count / f_bar

        return min(enc, ENC_MAX)

    @staticmethod
    def _expected_enc(s: float) -> float:
        # Wright (1990) expected ENC under mutation-only bias:
        # ENC = 2 + s + 29 / (s² + (1-s)²)
        denom = s ** 2 + (1 - s) ** 2
        if denom == 0:
            return ENC_MAX
        return min(2 + s + 29 / denom, ENC_MAX)

    def _make_plot(
        self,
        per_gene: list[dict[str, Any]],
        synonymous_families: dict[str, list[str]],
    ) -> None:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        gc3 = [row["GC3"] for row in per_gene if row["ENC"] is not None]
        enc = [row["ENC"] for row in per_gene if row["ENC"] is not None]

        xs = np.linspace(0.0, 1.0, 200)
        expected = [self._expected_enc(float(x)) for x in xs]

        fig, ax = plt.subplots()
        ax.scatter(gc3, enc, s=18, alpha=0.6, edgecolor="none", label="genes")
        ax.plot(xs, expected, color="black", linewidth=1.5, label="expected (no selection)")
        ax.set_xlabel("GC3")
        ax.set_ylabel("ENC")
        ax.set_title("ENC vs. GC3")
        ax.set_xlim(0, 1)
        ax.set_ylim(20, 65)
        ax.legend(loc="best")
        fig.tight_layout()

        out = self.output_file if self.output_file else "enc_plot.png"
        fig.savefig(out, dpi=150)
        plt.close(fig)

    def process_args(self, args: Any) -> dict[str, Any]:
        if args.translation_table is None:
            translation_table = "1"
        else:
            translation_table = args.translation_table
        return dict(
            fasta=args.fasta,
            translation_table=translation_table,
            verbose=bool(getattr(args, "verbose", False)),
            output_format=getattr(args, "format", None),
        )
