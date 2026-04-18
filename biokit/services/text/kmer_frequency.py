from collections import Counter
from typing import Any

from .base import Text
from ...helpers.files import read_and_parse_fasta_seqio


COMPLEMENT = str.maketrans("ACGT", "TGCA")


class KmerFrequency(Text):
    def __init__(self, args: Any) -> None:
        self.kmer_size: int = int(args.kmer_size)
        if self.kmer_size < 1:
            raise ValueError("kmer_size must be >= 1")
        self.canonical: bool = bool(getattr(args, "canonical", False))
        self.verbose: bool = bool(getattr(args, "verbose", False))
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        if self.fasta is None:
            raise ValueError("fasta cannot be None")
        output_format = self.normalize_output_format(self.output_format)
        records = read_and_parse_fasta_seqio(self.fasta)

        if self.verbose:
            rows: list[dict[str, Any]] = []
            for record in records:
                counts = self._count_kmers(str(record.seq))
                for kmer, count in self._sort_counts(counts):
                    rows.append({"id": record.id, "kmer": kmer, "count": count})
            field_order = ["id", "kmer", "count"]
            if output_format == "tsv":
                print("\t".join(field_order))
                for row in rows:
                    print("\t".join(str(row[f]) for f in field_order))
                return
            print(self.format_rows(rows, output_format, field_order=field_order))
            return

        aggregate: Counter[str] = Counter()
        for record in records:
            aggregate.update(self._count_kmers(str(record.seq)))

        field_order = ["kmer", "count"]
        if output_format == "tsv":
            print("\t".join(field_order))
            for kmer, count in self._sort_counts(aggregate):
                print(f"{kmer}\t{count}")
            return
        rows = [
            {"kmer": kmer, "count": count}
            for kmer, count in self._sort_counts(aggregate)
        ]
        print(self.format_rows(rows, output_format, field_order=field_order))

    def _count_kmers(self, sequence: str) -> Counter[str]:
        seq = sequence.upper()
        k = self.kmer_size
        counts: Counter[str] = Counter()
        valid = set("ACGT")
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            if not set(kmer).issubset(valid):
                continue
            if self.canonical:
                rc = kmer.translate(COMPLEMENT)[::-1]
                if rc < kmer:
                    kmer = rc
            counts[kmer] += 1
        return counts

    @staticmethod
    def _sort_counts(counts: Counter[str]) -> list[tuple[str, int]]:
        return sorted(counts.items(), key=lambda item: (-item[1], item[0]))

    def process_args(self, args: Any) -> dict[str, str | None]:
        return dict(fasta=args.fasta, output_format=getattr(args, "format", None))
