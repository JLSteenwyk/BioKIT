import sys
import textwrap
from typing import Any

from Bio import SeqIO

from .base import Text
from ...helpers.files import read_and_parse_fasta_seqio


class RenameFastaEntries(Text):
    def __init__(self, args: Any) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        if self.fasta is None or self.idmap is None or self.output_file_path is None:
            raise ValueError("fasta, idmap, and output_file_path cannot be None")
        # create biopython object of sequences
        records = read_and_parse_fasta_seqio(self.fasta)

        # save idmap to a dictionary
        idmap = self.idmap_to_dictionary(self.idmap)

        # replace and write out
        self.replace_ids_and_write(self.output_file_path, records, idmap)

    def process_args(self, args: Any) -> dict[str, str]:
        if args.output is None:
            output_file_path = f"{args.fasta}.renamed.fa"
        else:
            output_file_path = f"{args.output}"

        if output_file_path == '-.renamed.fa':
            output_file_path = 'renamed.fa'

        return dict(
            fasta=args.fasta, idmap=args.idmap, output_file_path=output_file_path
        )

    def replace_ids_and_write(self, output_file_path: str, records: Any, idmap: dict[str, str]) -> None:
        """
        for tips with a name as a key in the idmap,
        replace that tip name with the value in the
        idmap and write to output file
        """
        with open(output_file_path, "w") as output_handle:
            try:
                for record in records:
                    if record.id in idmap:
                        # replace ID
                        record.id = idmap[record.id]
                        # remove description
                        record.description = ""
                    SeqIO.write(record, output_handle, "fasta")
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

    def idmap_to_dictionary(self, idmap: str) -> dict[str, str]:
        """
        read idmap into a dictionary
        """
        id_map: dict[str, str] = {}
        try:
            with open(idmap) as identifiers:
                for line in identifiers:
                    stripped = line.strip()
                    if not stripped or stripped.startswith("#"):
                        continue
                    parts = stripped.split()
                    if len(parts) < 2:
                        raise ValueError(
                            f"Malformed idmap line in {idmap!r}: {line.rstrip()!r}"
                        )
                    key, val = parts[0], parts[1]
                    id_map[key] = val
            return id_map
        except FileNotFoundError:
            print(
                textwrap.dedent(
                    f"""
                    {idmap} corresponds to no such file or directory.
                    Please double check pathing and filenames
                    """
                )
            )
            sys.exit()
