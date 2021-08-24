from .base import Genome


class NumberOfLargeScaffolds(Genome):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        cnt, total_len = self.number_of_large_scaffolds()
        print(f"{cnt}\t{total_len}")

    def process_args(self, args):
        if args.threshold is None:
            threshold = 500
        else:
            threshold = int(args.threshold)

        return dict(
            fasta=args.fasta,
            threshold=threshold,
        )
