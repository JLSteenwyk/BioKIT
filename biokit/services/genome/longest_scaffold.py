from .base import Genome


class LongestScaffold(Genome):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        max_len = self.longest_scaffold()

        print(max_len)

    def process_args(self, args):
        return dict(fasta=args.fasta)
