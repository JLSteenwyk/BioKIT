import sys
from typing import Any

from Bio import SeqIO

from .base import Text


class FileFormatConverter(Text):
    def __init__(self, args: Any) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        file_formats = [
            "fasta",
            "clustal",
            "maf",
            "mauve",
            "phylip",
            "phylip_sequential",
            "phylip_relaxed",
            "stockholm",
        ]

        input_file_format = self.input_file_format
        output_file_format = self.output_file_format
        if input_file_format in file_formats and output_file_format in file_formats:
            input_file_format = self._normalize_format(input_file_format)
            output_file_format = self._normalize_format(output_file_format)

            SeqIO.convert(
                self.input_file,
                input_file_format,
                self.output_file,
                output_file_format,
            )
        else:
            print(
                f"File format not acceptable. Please use one of the following: {file_formats}"
            )

    @staticmethod
    def _normalize_format(file_format: str) -> str:
        if file_format == "phylip_sequential":
            return "phylip-sequential"
        if file_format == "phylip_relaxed":
            return "phylip-relaxed"
        return file_format

    def process_args(self, args: Any) -> dict[str, Any]:
        if args.input_file == '-':
            input_file = sys.stdin
        else:
            input_file = args.input_file
        return dict(
            input_file=input_file,
            output_file_format=args.output_file_format,
            input_file_format=args.input_file_format,
            output_file=args.output_file,
        )
