import sys
from typing import Any

from Bio.Seq import Seq

from .base import Text
from ...helpers.files import read_and_parse_fasta_seqio


class FindOrfs(Text):
    def __init__(self, args: Any) -> None:
        min_length = getattr(args, "min_length", None)
        self.min_length: int = int(min_length) if min_length is not None else 100
        if self.min_length < 1:
            raise ValueError("min_length must be >= 1")
        table = getattr(args, "translation_table", None)
        self.table: int = int(table) if table is not None else 1
        self.extract: bool = bool(getattr(args, "extract", False))
        self.protein: bool = bool(getattr(args, "protein", False))
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        if self.fasta is None:
            raise ValueError("fasta cannot be None")
        output_format = self.normalize_output_format(self.output_format)
        records = read_and_parse_fasta_seqio(self.fasta)

        orfs = []
        for record in records:
            orfs.extend(self._find_orfs_in_record(record.id, str(record.seq)))

        if self.extract:
            self._print_extracted_fasta(orfs)
            return

        field_order = ["id", "orf_id", "frame", "start", "stop", "length_nt", "length_aa"]
        rows = [
            {k: orf[k] for k in field_order}
            for orf in orfs
        ]
        if output_format == "tsv":
            print("\t".join(field_order))
            for row in rows:
                print("\t".join(str(row[f]) for f in field_order))
            return
        print(self.format_rows(rows, output_format, field_order=field_order))

    def _find_orfs_in_record(self, seq_id: str, sequence: str) -> list[dict[str, Any]]:
        seq = Seq(sequence.upper())
        seq_len = len(seq)
        results: list[dict[str, Any]] = []
        strand_pairs = [("+", seq), ("-", seq.reverse_complement())]

        for strand, nuc in strand_pairs:
            for frame in range(3):
                trimmed = nuc[frame:]
                usable_len = (len(trimmed) // 3) * 3
                if usable_len < 3:
                    continue
                trans = str(trimmed[:usable_len].translate(table=self.table))
                aa_pos = 0
                while aa_pos < len(trans):
                    next_stop = trans.find("*", aa_pos)
                    segment_end = next_stop if next_stop != -1 else len(trans)
                    m_pos = trans.find("M", aa_pos, segment_end)
                    if m_pos != -1:
                        aa_len = segment_end - m_pos
                        if aa_len >= self.min_length:
                            orf_nt = self._extract_orf_sequence(
                                nuc, frame, m_pos, segment_end, next_stop != -1
                            )
                            start, stop = self._compute_positions(
                                strand, seq_len, frame, m_pos, segment_end,
                                next_stop != -1,
                            )
                            results.append({
                                "id": seq_id,
                                "orf_id": f"{seq_id}_orf{len(results) + 1}",
                                "frame": f"{strand}{frame + 1}",
                                "start": start,
                                "stop": stop,
                                "length_nt": len(orf_nt),
                                "length_aa": aa_len,
                                "nt_seq": str(orf_nt),
                            })
                    aa_pos = segment_end + 1

        return results

    @staticmethod
    def _extract_orf_sequence(
        nuc: Seq, frame: int, m_pos: int, segment_end: int, has_stop: bool
    ) -> Seq:
        nt_start = frame + m_pos * 3
        nt_end = frame + segment_end * 3
        if has_stop:
            nt_end += 3
        return nuc[nt_start:nt_end]

    @staticmethod
    def _compute_positions(
        strand: str,
        seq_len: int,
        frame: int,
        m_pos: int,
        segment_end: int,
        has_stop: bool,
    ) -> tuple[int, int]:
        nt_start = frame + m_pos * 3
        nt_end = frame + segment_end * 3
        if has_stop:
            nt_end += 3
        if strand == "+":
            return nt_start + 1, nt_end
        return seq_len - nt_end + 1, seq_len - nt_start

    def _print_extracted_fasta(self, orfs: list[dict[str, Any]]) -> None:
        out = sys.stdout
        for orf in orfs:
            nt_seq = orf["nt_seq"]
            if self.protein:
                protein = str(
                    Seq(nt_seq).translate(table=self.table, to_stop=True)
                )
                seq_str = protein
            else:
                seq_str = nt_seq
            header = (
                f"{orf['orf_id']} source={orf['id']} frame={orf['frame']} "
                f"start={orf['start']} stop={orf['stop']} "
                f"length_nt={orf['length_nt']} length_aa={orf['length_aa']}"
            )
            out.write(f">{header}\n{seq_str}\n")

    def process_args(self, args: Any) -> dict[str, str | None]:
        return dict(fasta=args.fasta, output_format=getattr(args, "format", None))
