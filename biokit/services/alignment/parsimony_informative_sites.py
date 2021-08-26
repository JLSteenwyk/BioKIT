from .base import Alignment
from ...helpers.files import read_alignment_alignio


class ParsimonyInformativeSites(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment = read_alignment_alignio(self.fasta)

        # alignment length
        aln_len = alignment.get_alignment_length()

        pis, _, _ = self.determine_pis_vs_cs(alignment, aln_len)

        print(pis)

    def process_args(self, args):
        return dict(fasta=args.fasta)
