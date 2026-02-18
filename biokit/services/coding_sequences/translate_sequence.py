from typing import Any

from .base import CodingSequence
from ...helpers.files import iter_fasta_entries


class TranslateSequence(CodingSequence):
    def __init__(self, args: Any) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        # get the translation table as a dictionary
        translation_table = self.read_translation_table(self.translation_table)
        if self.fasta is None:
            raise ValueError("fasta cannot be None")
        if self.output_file_path is None:
            raise ValueError("output_file_path cannot be None")
        translation_table_get = translation_table.get

        with open(self.output_file_path, "w") as output_file_path:
            for seq_id, sequence in iter_fasta_entries(self.fasta):
                sequence = sequence.upper().replace("T", "U")
                amino_acids = []
                if len(sequence) % 3 == 0:
                    for position in range(0, len(sequence), 3):
                        codon = sequence[position:position + 3]
                        amino_acids.append(translation_table_get(codon, "X"))
                translated_sequence = "".join(amino_acids)
                if translated_sequence:
                    wrapped = "\n".join(
                        translated_sequence[i:i + 60]
                        for i in range(0, len(translated_sequence), 60)
                    )
                    output_file_path.write(f">{seq_id}\n{wrapped}\n")
                else:
                    output_file_path.write(f">{seq_id}\n\n")

    def process_args(self, args: Any) -> dict[str, str]:
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
