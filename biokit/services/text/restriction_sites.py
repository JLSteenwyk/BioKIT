from typing import Any

from Bio.Restriction import RestrictionBatch, Analysis  # type: ignore[attr-defined]
from Bio.Seq import Seq

from .base import Text
from ...helpers.files import read_and_parse_fasta_seqio


class RestrictionSites(Text):
    def __init__(self, args: Any) -> None:
        self.enzymes: list[str] = self.process_enzymes(args)
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        if self.fasta is None:
            raise ValueError("fasta cannot be None")

        rb = RestrictionBatch(self.enzymes)
        records = read_and_parse_fasta_seqio(self.fasta)

        for record in records:
            seq = Seq(str(record.seq))
            ana = Analysis(rb, seq)
            result = ana.full()

            for enzyme in sorted(result, key=lambda e: str(e)):
                positions = result[enzyme]
                if not positions:
                    print(f"{record.id}\t{enzyme}\tno_sites")
                    continue

                pos_str = ",".join(str(p) for p in positions)
                cuts = [0] + sorted(positions) + [len(seq)]
                fragments = [
                    cuts[i + 1] - cuts[i] for i in range(len(cuts) - 1)
                ]
                frag_str = ",".join(str(f) for f in fragments)
                print(
                    f"{record.id}\t{enzyme}\t{len(positions)}\t"
                    f"{pos_str}\t{frag_str}"
                )

    @staticmethod
    def process_enzymes(args: Any) -> list[str]:
        raw = args.enzymes
        enzymes: list[str] = []
        for item in raw:
            for e in item.split(","):
                e = e.strip()
                if e:
                    enzymes.append(e)
        return enzymes

    def process_args(self, args: Any) -> dict[str, str | None]:
        return dict(fasta=args.fasta)
