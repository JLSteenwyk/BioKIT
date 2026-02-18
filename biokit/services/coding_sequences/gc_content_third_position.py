from typing import Any

from .base import CodingSequence


class GCContentThirdPosition(CodingSequence):
    def __init__(self, args: Any) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        for line in self.gc_content_by_codon_position(2):
            print(line)

    def process_args(self, args: Any) -> dict[str, Any]:
        return dict(fasta=args.fasta, verbose=args.verbose)
