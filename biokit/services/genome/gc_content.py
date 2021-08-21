import re
import sys

from Bio import SeqIO

from .base import Genome


class GCContent(Genome):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        # create biopython object of sequences
        records = SeqIO.parse(self.fasta, "fasta")

        # initialize and populate dict for
        # holding the entry sequences
        entry_and_seq = {}
        for record in records:
            entry_and_seq[record.id] = str(record.seq)

        if self.verbose:
            for entry, seq in entry_and_seq.items():
                seq, matches = self.find_matches_and_remove_gaps(seq)
                try:
                    print(f"{entry}\t{round(len(matches)/len(seq), 4)}")
                except BrokenPipeError:
                    pass
        else:
            all_seqs = []
            for entry, seq in entry_and_seq.items():
                all_seqs.append(seq)
            all_seqs = "".join(all_seqs)
            all_seqs, matches = self.find_matches_and_remove_gaps(all_seqs)
            try:
                gc_content = round(len(matches) / len(all_seqs), 4)
            except ZeroDivisionError:
                try:
                    print(
                        "Input file has an unacceptable format. Please check input file argument."
                    )
                    sys.exit()
                except BrokenPipeError:
                    pass
            try:
                print(gc_content)
            except BrokenPipeError:
                pass

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
        )
