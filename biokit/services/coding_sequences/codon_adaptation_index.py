import math
from collections import Counter, defaultdict
from typing import Any

from .base import CodingSequence
from ...helpers.files import iter_fasta_entries


PSEUDOCOUNT = 0.5


class CodonAdaptationIndex(CodingSequence):
    def __init__(self, args: Any) -> None:
        self.reference: str = args.reference
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        if self.fasta is None:
            raise ValueError("fasta cannot be None")
        if not self.reference:
            raise ValueError("reference fasta is required")
        output_format = self.normalize_output_format(self.output_format)
        translation_table = self.read_translation_table(self.translation_table)
        synonymous_families = self._synonymous_families(translation_table)

        w_values = self._compute_w_values(synonymous_families)
        per_gene = self._compute_per_gene(translation_table, w_values)

        if output_format == "tsv":
            for row in per_gene:
                if self.verbose:
                    print(f"{row['gene_id']}\t{row['CAI']}\t{row['n_codons']}")
                else:
                    print(f"{row['gene_id']}\t{row['CAI']}")
            return

        if self.verbose:
            field_order = ["gene_id", "CAI", "n_codons"]
        else:
            field_order = ["gene_id", "CAI"]
        rows = [{k: row[k] for k in field_order} for row in per_gene]
        print(self.format_rows(rows, output_format, field_order=field_order))

    @staticmethod
    def _synonymous_families(
        translation_table: dict[str, str],
    ) -> dict[str, list[str]]:
        families: dict[str, list[str]] = defaultdict(list)
        for codon, aa in translation_table.items():
            if aa == "*":
                continue
            families[aa].append(codon)
        return dict(families)

    def _count_reference_codons(
        self,
        valid_codons: set[str],
    ) -> Counter[str]:
        counts: Counter[str] = Counter()
        for _, seq in iter_fasta_entries(self.reference):
            sequence = seq.upper().replace("T", "U")
            if len(sequence) % 3 != 0:
                continue
            for i in range(0, len(sequence), 3):
                codon = sequence[i:i + 3]
                if codon in valid_codons:
                    counts[codon] += 1
        return counts

    def _compute_w_values(
        self,
        synonymous_families: dict[str, list[str]],
    ) -> dict[str, float]:
        valid_codons = {
            codon
            for codons in synonymous_families.values()
            for codon in codons
        }
        reference_counts = self._count_reference_codons(valid_codons)
        w: dict[str, float] = {}
        for codons in synonymous_families.values():
            if len(codons) == 1:
                continue
            adjusted = {
                c: (reference_counts.get(c, 0) if reference_counts.get(c, 0) > 0
                    else PSEUDOCOUNT)
                for c in codons
            }
            max_count = max(adjusted.values())
            if max_count == 0:
                continue
            for c in codons:
                w[c] = adjusted[c] / max_count
        return w

    def _compute_per_gene(
        self,
        translation_table: dict[str, str],
        w_values: dict[str, float],
    ) -> list[dict[str, Any]]:
        assert self.fasta is not None
        valid_codons = set(translation_table.keys())
        results: list[dict[str, Any]] = []
        for seq_id, seq in iter_fasta_entries(self.fasta):
            sequence = seq.upper().replace("T", "U")
            if len(sequence) % 3 != 0:
                continue
            log_w_values: list[float] = []
            for i in range(0, len(sequence), 3):
                codon = sequence[i:i + 3]
                if codon not in valid_codons:
                    continue
                if translation_table[codon] == "*":
                    continue
                if codon not in w_values:
                    continue
                wi = w_values[codon]
                if wi <= 0:
                    continue
                log_w_values.append(math.log(wi))
            if not log_w_values:
                cai_value: float | None = None
            else:
                cai_value = round(
                    math.exp(sum(log_w_values) / len(log_w_values)), 4
                )
            results.append({
                "gene_id": seq_id,
                "CAI": cai_value,
                "n_codons": len(log_w_values),
            })
        return results

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
