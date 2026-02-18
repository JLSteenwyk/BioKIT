from typing import Any

from .base import Alignment
from ...helpers.files import read_alignment_alignio


class VariableSites(Alignment):
    def __init__(self, args: Any) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        if self.fasta is None:
            raise ValueError("fasta cannot be None")
        output_format = self.normalize_output_format(self.output_format)
        alignment = read_alignment_alignio(self.fasta)

        # alignment length
        aln_len = alignment.get_alignment_length()

        _, vs, _, site_summary = self.determine_pis_vs_cs(alignment, aln_len)

        if output_format == "tsv" and self.verbose:
            for i in site_summary:
                print(f"{i[0]}\t{i[1]}")
        elif output_format == "tsv":
            print(vs)
        elif self.verbose:
            rows = [{"site_index": i[0], "classification": i[1]} for i in site_summary]
            print(self.format_rows(rows, output_format))
        else:
            print(self.format_object({"variable_sites": vs}, output_format))

    def process_args(self, args: Any) -> dict[str, Any | None]:
        return dict(
            fasta=args.fasta,
            verbose=args.verbose,
            output_format=getattr(args, "format", None),
        )
