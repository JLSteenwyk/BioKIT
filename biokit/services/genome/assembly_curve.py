from typing import Any

from .base import Genome


class AssemblyCurve(Genome):
    def __init__(self, args: Any) -> None:
        self.output_file: str | None = getattr(args, "output", None)
        self.plot: bool = getattr(args, "plot", False)
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        if self.fasta is None:
            raise ValueError("fasta cannot be None")
        output_format = self.normalize_output_format(self.output_format)

        contig_lengths = sorted(self._get_contig_lengths(), reverse=True)
        cumulative = []
        running = 0
        for length in contig_lengths:
            running += length
            cumulative.append(running)

        ranks = list(range(1, len(contig_lengths) + 1))

        if self.plot:
            self._make_plot(ranks, cumulative)

        if output_format == "tsv":
            print("rank\tlength\tcumulative_length")
            for rank, length, cum in zip(ranks, contig_lengths, cumulative):
                print(f"{rank}\t{length}\t{cum}")
            return

        rows = [
            {"rank": rank, "length": length, "cumulative_length": cum}
            for rank, length, cum in zip(ranks, contig_lengths, cumulative)
        ]
        print(self.format_rows(rows, output_format))

    def _make_plot(self, ranks: list[int], cumulative: list[int]) -> None:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        ax.step(ranks, cumulative, where="post", linewidth=1.5)
        ax.set_xlabel("Contig rank")
        ax.set_ylabel("Cumulative length (bp)")
        ax.set_title("Assembly curve")
        fig.tight_layout()

        out = self.output_file if self.output_file else "assembly_curve.png"
        fig.savefig(out, dpi=150)
        plt.close(fig)

    def process_args(self, args: Any) -> dict[str, str | None]:
        return dict(
            fasta=args.fasta,
            output_format=getattr(args, "format", None),
        )
