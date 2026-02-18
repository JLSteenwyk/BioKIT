from typing import Any

from Bio import SeqIO

from .base import Text
from ...helpers.files import read_and_parse_fasta_seqio


class ReorderBySequenceLength(Text):
    def __init__(self, args: Any) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        if self.fasta is None or self.output_file_path is None:
            raise ValueError("fasta and output_file_path cannot be None")
        records = list(read_and_parse_fasta_seqio(self.fasta))
        records.sort(key=lambda record: len(record), reverse=True)
        with open(self.output_file_path, "w") as output_file_path:
            for record in records:
                SeqIO.write(record, output_file_path, "fasta")

    def process_args(self, args: Any) -> dict[str, str]:
        if args.output is None:
            output_file_path = f"{args.fasta}.reordered.fa"
        else:
            output_file_path = f"{args.output}"

        return dict(fasta=args.fasta, output_file_path=output_file_path)
