import random
import sys
from typing import Any

from Bio.SeqIO.FastaIO import SimpleFastaParser

from .base import Text


class SampleSequences(Text):
    def __init__(self, args: Any) -> None:
        self.number: int | None = (
            int(args.number) if getattr(args, "number", None) is not None else None
        )
        self.percent: float | None = (
            float(args.percent) if getattr(args, "percent", None) is not None else None
        )
        if self.number is not None and self.percent is not None:
            raise ValueError("specify either --number or --percent, not both")
        if self.number is not None and self.number < 0:
            raise ValueError("--number must be >= 0")
        if self.percent is not None and not 0 <= self.percent <= 100:
            raise ValueError("--percent must be between 0 and 100")
        self.seed: int | None = (
            int(args.seed) if getattr(args, "seed", None) is not None else None
        )
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        if self.fasta is None:
            raise ValueError("fasta cannot be None")

        entries = self._read_entries()
        total = len(entries)
        count = self._sample_count(total)

        rng = random.Random(self.seed)
        indices = sorted(rng.sample(range(total), count)) if count > 0 else []

        out_handle = (
            open(self.output_file, "w") if self.output_file else sys.stdout
        )
        try:
            for i in indices:
                title, seq = entries[i]
                out_handle.write(f">{title}\n{seq}\n")
        finally:
            if self.output_file:
                out_handle.close()

    def _read_entries(self) -> list[tuple[str, str]]:
        assert self.fasta is not None
        if self.fasta == "-":
            return list(SimpleFastaParser(sys.stdin))
        with open(self.fasta) as handle:
            return list(SimpleFastaParser(handle))

    def _sample_count(self, total: int) -> int:
        if self.number is not None:
            if self.number > total:
                raise ValueError(
                    f"requested {self.number} sequences but only {total} available"
                )
            return self.number
        percent = self.percent if self.percent is not None else 10.0
        return round(total * percent / 100)

    def process_args(self, args: Any) -> dict[str, Any]:
        return dict(
            fasta=args.fasta,
            output_file=getattr(args, "output", None),
        )
