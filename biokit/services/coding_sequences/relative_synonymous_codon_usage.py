from Bio import SeqIO

from .base import CodingSequence


class RelativeSynonymousCodonUsage(CodingSequence):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        translation_table = self.read_translation_table(self.translation_table)  # noqa

        # get codon_table
        codon_table = dict()
        for k, v in translation_table.items():
            if v in codon_table.keys():
                codon_table[v].append(k)
            else:
                codon_table[v] = [k]

        # get counts of codons
        codon_counts = dict()
        records = SeqIO.parse(self.fasta, "fasta")
        for seq_record in records:
            if len(seq_record._seq) % 3 == 0:
                for position in range(0, len(seq_record._seq), 3):
                    codon = (
                        seq_record._seq[position : position + 3]
                        ._data.upper()
                        .replace("T", "U")
                    )
                    if codon in codon_counts.keys():
                        codon_counts[codon] += 1
                    else:
                        codon_counts[codon] = 1

        # calculate rscu
        rscu = dict()
        for _, codons in codon_table.items():
            observed_sum = 0
            for codon in codons:
                try:
                    observed_sum += codon_counts[codon]
                except KeyError:
                    observed_sum = 0
            for codon in codons:
                try:
                    rscu[codon] = round(
                        codon_counts[codon] / (observed_sum / len(codons)), 4
                    )
                except ZeroDivisionError:
                    rscu[codon] = 0
                except KeyError:
                    rscu[codon] = 0

        # reverse sort according to rscu values
        rscu = dict(sorted(rscu.items(), key=lambda item: item[1], reverse=True))

        res = ""
        i = 1
        num_codons = len(rscu)
        for codon, rscu in rscu.items():
            if i != num_codons:
                res += f"{codon}\t{rscu}\n"
                i += 1
            else:
                res += f"{codon}\t{rscu}"

        print(res)

    def process_args(self, args):
        if args.translation_table is None:
            translation_table = "1"
        else:
            translation_table = args.translation_table

        return dict(fasta=args.fasta, translation_table=translation_table)
