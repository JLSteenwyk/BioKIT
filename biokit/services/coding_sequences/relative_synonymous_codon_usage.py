from typing import Any

from .base import CodingSequence


class RelativeSynonymousCodonUsage(CodingSequence):
    def __init__(self, args: Any) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        output_format = self.normalize_output_format(self.output_format)
        translation_table = self.read_translation_table(self.translation_table)  # noqa

        # get rscu values
        rscu = self.calculate_rscu(translation_table)

        # reverse sort according to rscu values (preserve historical
        # tie ordering for default TSV output).
        sorted_items = sorted(rscu.items(), key=lambda item: item[1], reverse=True)

        if output_format == "tsv":
            print("\n".join(f"{codon}\t{rscu_value}" for codon, rscu_value in sorted_items))
            return

        deterministic_rows = sorted(
            sorted_items, key=lambda item: (-float(item[1]), item[0])
        )
        rows = [{"codon": codon, "rscu": rscu_value} for codon, rscu_value in deterministic_rows]
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
