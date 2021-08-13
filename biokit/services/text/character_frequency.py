import sys

from .base import Text

class CharacterFrequency(Text):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        char_freq = self.character_frequency()
        sum_of_chars = sum(char_freq.values())

        for char, freq in char_freq.items():
            print(f"{char}\t{round(freq/sum_of_chars, 4)}")
    
    def process_args(self, args):
        return dict(fasta=args.fasta)