import sys

from ..base import BaseService
from ...helpers.files import get_alignment_and_format as get_alignment_and_format_helper

class Text(BaseService):
    def __init__(
        self,
        *args,
        alignment_file_path=None,
        fasta=None,
        output_file_path=None,
        protein_file_path=None,
        nucleotide_file_path=None,
        alignment_list_path=None,
        prefix=None,
        idmap=None,
        reference=None,
        verbose=None,
        entry=None,
        input_file=None,
        output_file_format=None,
        input_file_format=None,
        output_file=None,
        threshold=None,
        ambiguous_character = None,
    ):
        self.alignment_file_path = alignment_file_path
        self.output_file_path = output_file_path
        self.protein_file_path = protein_file_path,
        self.nucleotide_file_path = nucleotide_file_path 
        self.alignment_list_path = alignment_list_path
        self.prefix = prefix
        self.fasta = fasta
        self.idmap = idmap
        self.reference = reference
        self.verbose = verbose
        self.entry = entry
        self.input_file = input_file
        self.output_file_format = output_file_format
        self.input_file_format = input_file_format
        self.output_file = output_file
        self.threshold = threshold
        self.ambiguous_character = ambiguous_character

    def get_alignment_and_format(self):
        """
        automatic file type determination
        """
        try:
            return get_alignment_and_format_helper(self.alignment_file_path)
        except FileNotFoundError:
            print("Input corresponds to no such file or directory.")
            print("Please double check pathing and filenames")
            sys.exit()

