from typing import Any

from Bio import SeqIO

from .base import Text


class Faidx(Text):
    def __init__(self, args: Any) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        if self.fasta is None or self.entry is None:
            raise ValueError("fasta and entry cannot be None")
        record_dict = SeqIO.index(self.fasta, "fasta")
        try:
            if self.entry not in record_dict:
                raise Exception(f"Entry {self.entry!r} not found in {self.fasta}")
            print(f">{record_dict[self.entry].name}\n{record_dict[self.entry].seq}")
        finally:
            record_dict.close()

    def process_args(self, args: Any) -> dict[str, str]:
        return dict(fasta=args.fasta, entry=args.entry)
