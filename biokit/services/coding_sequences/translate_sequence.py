from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .base import CodingSequence
from ...helpers.files import read_and_parse_fasta_seqio


class TranslateSequence(CodingSequence):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        # get the translation table as a dictionary
        translation_table = self.read_translation_table(self.translation_table)
        records = read_and_parse_fasta_seqio(self.fasta)

        with open(self.output_file_path, "w") as output_file_path:
            for seq_record in records:
                amino_acids = []
                if len(seq_record._seq) % 3 == 0:
                    for position in range(0, len(seq_record._seq), 3):
                        codon = (
                            seq_record._seq[position:position + 3]
                            ._data
                            .upper()
                            .replace(b"T", b"U")
                            .decode("utf-8")
                        )
                        if codon in translation_table.keys():
                            amino_acids.append(translation_table[codon])
                        else:
                            amino_acids.append("X")
                translated_seq_record = SeqRecord(
                    Seq("".join(amino_acids)),
                    id=seq_record.id,
                    name="",
                    description="",
                )

                SeqIO.write(translated_seq_record, output_file_path, "fasta")

    def process_args(self, args):
        if args.output is None:
            output_file_path = f"{args.fasta}.translated.fa"            
        else:
            output_file_path = f"{args.output}"

        if output_file_path == "-.translated.fa":
            output_file_path = "translated_seq.fa"

        if args.translation_table is None:
            translation_table = "1"
        else:
            translation_table = args.translation_table

        return dict(
            fasta=args.fasta,
            translation_table=translation_table,
            output_file_path=output_file_path,
        )
