from typing import Any

from .base import Text
from ...helpers.calculate_character_frequencies import calculate_character_frequencies


class CharacterFrequency(Text):
    def __init__(self, args: Any) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        if self.fasta is None:
            raise ValueError("fasta cannot be None")
        output_format = self.normalize_output_format(self.output_format)
        char_freq = calculate_character_frequencies(self.fasta)
        sum_of_chars = sum(char_freq.values())

        if output_format == "tsv":
            for char, freq in char_freq.items():
                print(f"{char}\t{round(freq/sum_of_chars, 4)}")
            return

        rows = [
            {"character": char, "frequency": round(freq / sum_of_chars, 4)}
            for char, freq in sorted(char_freq.items(), key=lambda item: item[0])
        ]
        print(self.format_rows(rows, output_format))

    def process_args(self, args: Any) -> dict[str, str | None]:
        return dict(fasta=args.fasta, output_format=getattr(args, "format", None))
