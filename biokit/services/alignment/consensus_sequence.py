import re
from collections import Counter

from Bio.Seq import Seq

from .base import Alignment
from ...helpers.files import read_alignment_alignio


def _dumb_consensus(alignment, threshold=0.7, ambiguous="X", require_multiple=False):
    consensus = ""
    con_len = alignment.get_alignment_length()

    for n in range(con_len):
        atom_dict = Counter()
        num_atoms = 0

        for record in alignment:
            try:
                residue = record[n]
            except IndexError:
                continue
            if residue != "-" and residue != ".":
                atom_dict[residue] += 1
                num_atoms += 1

        max_atoms = []
        max_size = 0
        for atom in atom_dict:
            if atom_dict[atom] > max_size:
                max_atoms = [atom]
                max_size = atom_dict[atom]
            elif atom_dict[atom] == max_size:
                max_atoms.append(atom)

        if require_multiple and num_atoms == 1:
            consensus += ambiguous
        elif len(max_atoms) == 1 and max_size / num_atoms >= threshold:
            consensus += max_atoms[0]
        else:
            consensus += ambiguous

    return Seq(consensus)


class ConsensusSequence(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment = read_alignment_alignio(self.fasta)

        if not self.ambiguous_character:
            ambiguous_character = "N"
        else:
            ambiguous_character = self.ambiguous_character

        if not self.threshold:
            threshold = 0.7
        else:
            threshold = float(self.threshold)

        consensus = _dumb_consensus(
            alignment,
            threshold=threshold, ambiguous=ambiguous_character
        )

        header = ">" + re.sub("^.*/", "", str(self.fasta)) + ".consensus"
        print(f"{header}\n{consensus}")

    def process_args(self, args):
        return dict(
            fasta=args.fasta,
            threshold=args.threshold,
            ambiguous_character=args.ambiguous_character,
        )
