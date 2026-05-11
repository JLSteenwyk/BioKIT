from typing import Any

import numpy as np

from .base import CodingSequence
from ...helpers.files import iter_fasta_entries


class NeutralityPlot(CodingSequence):
    def __init__(self, args: Any) -> None:
        self.plot: bool = bool(getattr(args, "plot", False))
        self.output_file: str | None = getattr(args, "output", None)
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        if self.fasta is None:
            raise ValueError("fasta cannot be None")
        output_format = self.normalize_output_format(self.output_format)

        per_gene = self._per_gene_gc12_gc3()
        regression = self._regression(per_gene)

        if self.plot:
            self._make_plot(per_gene, regression)

        if output_format == "tsv":
            if self.verbose:
                print("id\tGC12\tGC3")
                for row in per_gene:
                    print(f"{row['id']}\t{row['GC12']}\t{row['GC3']}")
                print()
            for label, value in regression.items():
                print(f"{value}\t{label}")
            return

        payload: dict[str, Any] = {"regression": regression}
        if self.verbose:
            payload["per_gene"] = per_gene
        print(self.format_object(payload, output_format))

    def _per_gene_gc12_gc3(self) -> list[dict[str, Any]]:
        assert self.fasta is not None
        rows: list[dict[str, Any]] = []
        for seq_id, seq in iter_fasta_entries(self.fasta):
            pos1 = self.get_codon_position_chars(seq, 0)
            pos2 = self.get_codon_position_chars(seq, 1)
            pos3 = self.get_codon_position_chars(seq, 2)
            gc1 = self.calculate_gc_content(pos1)
            gc2 = self.calculate_gc_content(pos2)
            gc3 = self.calculate_gc_content(pos3)
            gc12 = round((gc1 + gc2) / 2, 4)
            rows.append({"id": seq_id, "GC12": gc12, "GC3": gc3})
        return rows

    @staticmethod
    def _regression(per_gene: list[dict[str, Any]]) -> dict[str, Any]:
        n = len(per_gene)
        if n < 2:
            return {
                "slope": None,
                "intercept": None,
                "r_squared": None,
                "n": n,
            }
        gc3 = np.array([row["GC3"] for row in per_gene], dtype=float)
        gc12 = np.array([row["GC12"] for row in per_gene], dtype=float)
        if gc3.std() == 0:
            return {
                "slope": 0.0,
                "intercept": round(float(gc12.mean()), 4),
                "r_squared": 0.0,
                "n": n,
            }
        slope, intercept = np.polyfit(gc3, gc12, 1)
        predicted = slope * gc3 + intercept
        ss_res = float(np.sum((gc12 - predicted) ** 2))
        ss_tot = float(np.sum((gc12 - gc12.mean()) ** 2))
        r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
        return {
            "slope": round(float(slope), 4),
            "intercept": round(float(intercept), 4),
            "r_squared": round(float(r_squared), 4),
            "n": n,
        }

    def _make_plot(
        self,
        per_gene: list[dict[str, Any]],
        regression: dict[str, Any],
    ) -> None:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        gc3 = [row["GC3"] for row in per_gene]
        gc12 = [row["GC12"] for row in per_gene]

        fig, ax = plt.subplots()
        ax.scatter(gc3, gc12, s=20, alpha=0.6, edgecolor="none")
        if regression["slope"] is not None and len(per_gene) >= 2:
            xs = np.array([min(gc3), max(gc3)])
            ys = regression["slope"] * xs + regression["intercept"]
            label = (
                f"slope={regression['slope']:.3f}, "
                f"R²={regression['r_squared']:.3f}"
            )
            ax.plot(xs, ys, color="black", linewidth=1.5, label=label)
            ax.legend(loc="best")
        ax.set_xlabel("GC3")
        ax.set_ylabel("GC12")
        ax.set_title("Neutrality plot")
        fig.tight_layout()

        out = self.output_file if self.output_file else "neutrality_plot.png"
        fig.savefig(out, dpi=150)
        plt.close(fig)

    def process_args(self, args: Any) -> dict[str, Any]:
        return dict(
            fasta=args.fasta,
            verbose=bool(getattr(args, "verbose", False)),
            output_format=getattr(args, "format", None),
        )
