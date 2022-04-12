from .base import Alignment
from ...helpers.files import read_alignment_alignio


class ConstantSites(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment = read_alignment_alignio(self.fasta)

        # alignment length
        aln_len = alignment.get_alignment_length()

        _, _, cs, site_summary = self.determine_pis_vs_cs(alignment, aln_len)

        if self.verbose:
            for i in site_summary:
                print(f"{i[0]}\t{i[1]}")
        else:
            print(cs)

    def process_args(self, args):
        return dict(fasta=args.fasta, verbose=args.verbose)
