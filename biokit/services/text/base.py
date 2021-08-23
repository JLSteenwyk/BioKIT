from os import path

from ..base import BaseService

here = path.dirname(__file__)


class Text(BaseService):
    def __init__(
        self,
        *args,
        ambiguous_character=None,
        entry=None,
        fasta=None,
        idmap=None,
        input_file=None,
        input_file_format=None,
        output=None,
        output_file=None,
        output_file_format=None,
        output_file_path=None,
        reverse=None,
        threshold=None,
    ):
        self.ambiguous_character = ambiguous_character
        self.entry = entry
        self.fasta = fasta
        self.idmap = idmap
        self.input_file = input_file
        self.input_file_format = input_file_format
        self.output = output
        self.output_file = output_file
        self.output_file_format = output_file_format
        self.output_file_path = output_file_path
        self.reverse = reverse
        self.threshold = threshold
