import sys
from typing import Any

from Bio.SeqIO.QualityIO import FastqGeneralIterator

from .base import FastQ


class FastQToFasta(FastQ):
    def __init__(self, args: Any) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        if self.fastq is None:
            raise ValueError("fastq cannot be None")

        out_handle = (
            open(self.output_file, "w") if self.output_file else sys.stdout
        )
        try:
            if self.fastq == "-":
                for title, seq, _ in FastqGeneralIterator(sys.stdin):
                    out_handle.write(f">{title}\n{seq}\n")
            else:
                with open(self.fastq) as in_handle:
                    for title, seq, _ in FastqGeneralIterator(in_handle):
                        out_handle.write(f">{title}\n{seq}\n")
        finally:
            if self.output_file:
                out_handle.close()

    def process_args(self, args: Any) -> dict[str, Any]:
        return dict(
            fastq=args.fastq,
            output_file=getattr(args, "output", None),
        )
