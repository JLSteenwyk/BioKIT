from typing import Any

from .base import CodingSequence


class GCContentFourFoldDegenerateSites(CodingSequence):
    def __init__(self, args: Any) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        trans_table = self.read_translation_table(self.translation_table)
        for line in self.gc_content_four_fold_degenerate(trans_table):
            print(line)

    def process_args(self, args: Any) -> dict[str, Any]:
        return dict(
            fasta=args.fasta,
            verbose=args.verbose,
            translation_table=getattr(args, "translation_table", None),
        )
