from .base import CodingSequence


class RelativeSynonymousCodonUsage(CodingSequence):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        translation_table = self.read_translation_table(self.translation_table)  # noqa

        # get rscu values
        rscu = self.calculate_rscu(translation_table)

        # reverse sort according to rscu values
        rscu = dict(sorted(rscu.items(), key=lambda item: item[1], reverse=True))

        res = ""
        i = 1
        num_codons = len(rscu)
        for codon, rscu in rscu.items():
            if i != num_codons:
                res += f"{codon}\t{rscu}\n"
                i += 1
            else:
                res += f"{codon}\t{rscu}"

        print(res)

    def process_args(self, args):
        if args.translation_table is None:
            translation_table = "1"
        else:
            translation_table = args.translation_table

        return dict(fasta=args.fasta, translation_table=translation_table)
