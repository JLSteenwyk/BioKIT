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

        input_file_format = self.input_file_format
        output_file_format = self.output_file_format
        if input_file_format and output_file_format in file_formats:
            if input_file_format == "phylip_sequential":
                input_file_format = "phylip-sequential"
            elif input_file_format == "phylip_relaxed":
                input_file_format = "phylip-relaxed"

            if output_file_format == "phylip_sequential":
                output_file_format = "phylip-sequential"
            elif output_file_format == "phylip_relaxed":
                output_file_format = "phylip-relaxed"

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

    def process_args(self, args):
        return dict(
            input_file=args.input_file,
            output_file_format=args.output_file_format,
            input_file_format=args.input_file_format,
            output_file=args.output_file,
        )
