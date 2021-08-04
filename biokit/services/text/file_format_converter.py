from Bio import SeqIO

from .base import Text

class FileFormatConverter(Text):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
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

        if self.input_file_format and self.output_file_format in file_formats:
            count = SeqIO.convert(
                self.input_file,
                self.input_file_format,
                self.output_file,
                self.output_file_format
            )
        else:
            print(f"File format not acceptable. Please use one of the following: {file_formats}")

    def process_args(self, args):
        return dict(
            input_file=args.input_file,
            output_file_format=args.output_file_format,
            input_file_format=args.input_file_format,
            output_file=args.output_file
        )