from typing import Any

from ..base import BaseService


class Alignment(BaseService):
    def __init__(self, *args: Any, **kwargs: Any) -> None:
        self.code: str | None = kwargs.get("code")
        self.fasta: str | None = kwargs.get("fasta")
        self.ambiguous_character: str | None = kwargs.get("ambiguous_character")
        self.threshold: str | float | None = kwargs.get("threshold")
        self.verbose: bool | None = kwargs.get("verbose")
        self.output_format: str | None = kwargs.get("output_format")

    def determine_pis_vs_cs(self, alignment: Any, alignment_length: int) -> tuple[int, int, int, list[list[Any]]]:
        """
        determine number of parsimony informative,
        variable, and constant sites in an alignment
        """
        parsimony_informative_sites = 0
        variable_sites = 0
        constant_sites = 0
        site_summary: list[list[Any]] = []

        for i in range(0, alignment_length):
            temp: list[Any] = []
            temp.append(i)
            seq_at_position = ""
            seq_at_position += alignment[:, i]
            seq_at_position = seq_at_position.upper()
            seq_at_position = seq_at_position.replace("?", "")
            seq_at_position = seq_at_position.replace("-", "")
            num_occurences = {}
            for char in set(seq_at_position):
                num_occurences[char] = seq_at_position.count(char)
            d = dict((k, v) for k, v in num_occurences.items() if v >= 2)

            # two characters that occur at least twice
            if len(d) >= 2:
                parsimony_informative_sites += 1
                temp.append("parsimony_informative_site")

            # if one character occurs at least twice and is the only character,
            # the site is not parismony informative but it is constant
            elif len(d) == 1 and len(num_occurences) >= 1:
                constant_sites += 1
                temp.append("constant_site")
            else:
                temp.append("Not_pis_vs_cs")

            if len(d) > 1 and len(num_occurences) >= 2:
                variable_sites += 1
                if temp[1]:
                    new_str = temp[1] + "_and_variable_site"
                    temp[1] = new_str
                else:
                    temp.append("variable_site")

            site_summary.append(temp)

        return parsimony_informative_sites, variable_sites, constant_sites, site_summary
