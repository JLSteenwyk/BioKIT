import statistics as stat

from Bio import SeqIO

from .base import CodingSequence


class GeneWiseRelativeSynonymousCodonUsage(CodingSequence):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        translation_table = self.read_translation_table(self.translation_table)  # noqa

        # get rscu values
        rscu = self.calculate_rscu(translation_table)

        records = SeqIO.parse(self.fasta, "fasta")

        gw_rscu = []
        for seq_record in records:
            if len(seq_record._seq) % 3 == 0:
                rscus_curr_gene = []
                for position in range(0, len(seq_record._seq), 3):
                    codon = (
                        seq_record._seq[position:position + 3]
                        ._data.upper()
                        .replace("T", "U")
                    )
                    rscus_curr_gene.append(rscu[codon])
                temp = []
                temp.append(seq_record.id)
                temp.append(round(stat.mean(rscus_curr_gene), 4))
                temp.append(round(stat.median(rscus_curr_gene), 4))
                temp.append(round(stat.stdev(rscus_curr_gene), 4))
                gw_rscu.append(temp)

        res = ""
        i = 1
        num_genes = len(gw_rscu)
        for gene_stats in gw_rscu:
            if i != num_genes:
                res += f"{gene_stats[0]}\t{gene_stats[1]}\t{gene_stats[2]}\t{gene_stats[3]}\n"
                i += 1
            else:
                res += f"{gene_stats[0]}\t{gene_stats[1]}\t{gene_stats[2]}\t{gene_stats[3]}"

        print(res)

    def process_args(self, args):
        if args.translation_table is None:
            translation_table = "1"
        else:
            translation_table = args.translation_table

        return dict(fasta=args.fasta, translation_table=translation_table)
