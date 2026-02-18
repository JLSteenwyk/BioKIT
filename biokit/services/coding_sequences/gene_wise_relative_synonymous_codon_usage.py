import statistics as stat
from typing import Any

from .base import CodingSequence
from ...helpers.files import iter_fasta_entries


class GeneWiseRelativeSynonymousCodonUsage(CodingSequence):
    def __init__(self, args: Any) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        output_format = self.normalize_output_format(self.output_format)
        translation_table = self.read_translation_table(self.translation_table)  # noqa

        # get rscu values
        rscu = self.calculate_rscu(translation_table)

        gw_rscu: list[tuple[str, float, float, float]] = []
        translation_table_keys = set(translation_table.keys())
        if self.fasta is None:
            raise ValueError("fasta cannot be None")
        for seq_id, seq in iter_fasta_entries(self.fasta):
            sequence = seq.upper().replace("T", "U")
            if len(sequence) % 3 != 0:
                continue
            rscus_curr_gene: list[float] = []
            for position in range(0, len(sequence), 3):
                codon = sequence[position:position + 3]
                if codon in translation_table_keys:
                    rscus_curr_gene.append(float(rscu[codon]))
            if not rscus_curr_gene:
                continue
            if len(rscus_curr_gene) >= 2:
                std_dev = round(stat.stdev(rscus_curr_gene), 4)
            else:
                std_dev = 0.0
            gw_rscu.append(
                (
                    seq_id,
                    round(stat.mean(rscus_curr_gene), 4),
                    round(stat.median(rscus_curr_gene), 4),
                    std_dev,
                )
            )

        if output_format == "tsv":
            print(
                "\n".join(
                    f"{gene_stats[0]}\t{gene_stats[1]}\t{gene_stats[2]}\t{gene_stats[3]}"
                    for gene_stats in gw_rscu
                )
            )
            return

        rows = [
            {
                "gene_id": gene_id,
                "mean_rscu": mean_rscu,
                "median_rscu": median_rscu,
                "stddev_rscu": stddev_rscu,
            }
            for gene_id, mean_rscu, median_rscu, stddev_rscu in gw_rscu
        ]
        rows.sort(key=lambda row: row["gene_id"])
        print(self.format_rows(rows, output_format))

    def process_args(self, args: Any) -> dict[str, str | None]:
        if args.translation_table is None:
            translation_table = "1"
        else:
            translation_table = args.translation_table

        return dict(
            fasta=args.fasta,
            translation_table=translation_table,
            output_format=getattr(args, "format", None),
        )
