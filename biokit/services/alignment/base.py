from ..base import BaseService


class Alignment(BaseService):
    def __init__(
        self,
        *args,
        code=None,
        fasta=None,
        ambiguous_character=None,
        threshold=None,
    ):
        self.code = code
        self.fasta = fasta
        self.ambiguous_character = ambiguous_character
        self.threshold = threshold

    def determine_pis_vs_cs(self, alignment, alignment_length):
        """
        determine number of parsimony informative,
        variable, and constant sites in an alignment
        """
        parsimony_informative_sites = 0
        variable_sites = 0
        constant_sites = 0
        for i in range(0, alignment_length):
            seq_at_position = ""
            seq_at_position += alignment[:, i]
            seq_at_position = seq_at_position.replace("?", "")
            seq_at_position = seq_at_position.replace("-", "")
            num_occurences = {}
            for char in set(seq_at_position):
                num_occurences[char] = seq_at_position.count(char)
            d = dict((k, v) for k, v in num_occurences.items() if v >= 2)

            # two characters that occur at least twice
            if len(d) >= 2:
                parsimony_informative_sites += 1
            # if one character occurs at least twice and is the only character,
            # the site is not parismony informative but it is constant
            elif len(d) == 1 and len(num_occurences) == 1:
                constant_sites += 1
            else:
                pass

            if len(num_occurences) >= 2:
                variable_sites += 1

        return parsimony_informative_sites, variable_sites, constant_sites
