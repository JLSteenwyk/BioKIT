import sys
import textwrap

from Bio import SeqIO

from .base import Text
from ...helpers.files import read_and_parse_fasta_seqio


class RenameFastaEntries(Text):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        # create biopython object of sequences
        records = read_and_parse_fasta_seqio(self.fasta)

        # save idmap to a dictionary
        idmap = self.idmap_to_dictionary(self.idmap)

        # replace and write out
        self.replace_ids_and_write(self.output_file_path, records, idmap)

    def process_args(self, args):
        if args.output is None:
            output_file_path = f"{args.fasta}.renamed.fa"
        else:
            output_file_path = f"{args.output}"

        if output_file_path == '-.renamed.fa':
            output_file_path = 'renamed.fa'

        return dict(
            fasta=args.fasta, idmap=args.idmap, output_file_path=output_file_path
        )

    def replace_ids_and_write(self, output_file_path, records, idmap):
        """
        for tips with a name as a key in the idmap,
        replace that tip name with the value in the
        idmap and write to output file
        """
        with open(output_file_path, "w") as output_file_path:
            try:
                for record in records:
                    if record.id in idmap:
                        # replace ID
                        record.id = idmap[record.id]
                        # remove description
                        record.description = ""
                    SeqIO.write(record, output_file_path, "fasta")
            except FileNotFoundError:
                print(
                    textwrap.dedent(
                        f"""
                        {self.fasta} corresponds to no such file or directory.
                        Please double check pathing and filenames
                        """
                    )
                )
                sys.exit()

    def idmap_to_dictionary(self, idmap: str) -> dict:
        """
        read idmap into a dictionary
        """
        idmap = {}
        try:
            with open(self.idmap) as identifiers:
                for line in identifiers:
                    (key, val) = line.split()
                    idmap[key] = val
            return idmap
        except FileNotFoundError:
            print(
                textwrap.dedent(
                    f"""
                    {self.idmap} corresponds to no such file or directory.
                    Please double check pathing and filenames
                    """
                )
            )
            sys.exit()
