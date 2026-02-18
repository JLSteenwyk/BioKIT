from os import path
from typing import Any

from ..base import BaseService

here = path.dirname(__file__)


class Text(BaseService):
    def __init__(self, *args: Any, **kwargs: Any) -> None:
        self.ambiguous_character: str | None = kwargs.get("ambiguous_character")
        self.entry: str | None = kwargs.get("entry")
        self.fasta: str | None = kwargs.get("fasta")
        self.idmap: str | None = kwargs.get("idmap")
        self.input_file: str | Any | None = kwargs.get("input_file")
        self.input_file_format: str | None = kwargs.get("input_file_format")
        self.output: str | None = kwargs.get("output")
        self.output_file: str | None = kwargs.get("output_file")
        self.output_file_format: str | None = kwargs.get("output_file_format")
        self.output_file_path: str | None = kwargs.get("output_file_path")
        self.output_format: str | None = kwargs.get("output_format")
        self.reverse: bool | None = kwargs.get("reverse")
        self.threshold: int | None = kwargs.get("threshold")
