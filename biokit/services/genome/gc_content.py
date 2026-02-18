import re

from .base import Genome
from ...helpers.files import read_and_parse_fasta_seqio


class GCContent(Genome):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        output_format = self.normalize_output_format(self.output_format)
        # create biopython object of sequences
        records = read_and_parse_fasta_seqio(self.fasta)

        # initialize and populate dict for
        # holding the entry sequences
        entry_and_seq = {}
        for record in records:
            entry_and_seq[record.id] = str(record.seq)

        if self.verbose:
            if output_format != "tsv":
                rows = []
            for entry, seq in entry_and_seq.items():
                seq, matches = self.find_matches_and_remove_gaps(seq)
                if len(seq) == 0:
                    gc_content = 0
                else:
                    gc_content = round(len(matches) / len(seq), 4)
                if output_format == "tsv":
                    print(f"{entry}\t{gc_content}")
                else:
                    rows.append({"entry": entry, "gc_content": gc_content})
            if output_format != "tsv":
                rows.sort(key=lambda row: row["entry"])
                print(self.format_rows(rows, output_format))
        else:
            all_seqs = []
            for entry, seq in entry_and_seq.items():
                all_seqs.append(seq)
            all_seqs = "".join(all_seqs)
            all_seqs, matches = self.find_matches_and_remove_gaps(all_seqs)
            if len(all_seqs) == 0:
                gc_content = 0
            else:
                gc_content = round(len(matches) / len(all_seqs), 4)
            if output_format == "tsv":
                print(gc_content)
            else:
                print(self.format_object({"gc_content": gc_content}, output_format))

    def find_matches_and_remove_gaps(self, seq: str):
        regex_pattern = re.compile("[GgCc]")
        seq = seq.replace("-", "")
        seq = seq.replace("?", "")
        matches = regex_pattern.findall(seq)
        return seq, matches

    def process_args(self, args):
        return dict(
            fasta=args.fasta,
            verbose=args.verbose,
            output_format=getattr(args, "format", None),
        )
