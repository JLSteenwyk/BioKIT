from ..base import BaseService


class CodingSequence(BaseService):
    def __init__(
        self,
        *args,
        fasta=None,
        verbose=None,
        translation_table=None,
        output_file_path=None,
    ):
        self.fasta = fasta
        self.verbose = verbose
        self.translation_table = translation_table
        self.output_file_path = output_file_path
