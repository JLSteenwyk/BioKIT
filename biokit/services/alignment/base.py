from ..base import BaseService


class Alignment(BaseService):
    def __init__(
        self,
        *args,
        code=None,
        fasta=None,
        ambiguous_character=None,
        threshold=None,
        verbose=None,
    ):
        self.code = code
        self.fasta = fasta
        self.ambiguous_character = ambiguous_character
        self.threshold = threshold
        self.verbose = verbose

    def determine_pis_vs_cs(self, alignment, alignment_length):
        """
        determine number of parsimony informative,
        variable, and constant sites in an alignment
        """
        parsimony_informative_sites = 0
        variable_sites = 0
        constant_sites = 0
        site_summary = []

        for i in range(0, alignment_length):
            temp = []
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
