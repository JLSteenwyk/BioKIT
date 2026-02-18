from collections import Counter
from typing import Any
from Bio.Seq import Seq

from .base import Alignment
from ...helpers.files import read_alignment_alignio


class PSSM:
    def __init__(self, pssm: list[tuple[str, dict[str, float]]]) -> None:
        self.pssm = pssm

    def __getitem__(self, pos: int) -> dict[str, float]:
        return self.pssm[pos][1]

    def __str__(self) -> str:
        out = " "
        all_residues = sorted(self.pssm[0][1])

        for res in all_residues:
            out += "   %s" % res
        out += "\n"

        for item in self.pssm:
            out += "%s " % item[0]
            for res in all_residues:
                out += " %.1f" % item[1][res]

            out += "\n"
        return out


def _dumb_consensus(
    alignment: Any,
    threshold: float = 0.7,
    ambiguous: str = "X",
    require_multiple: bool = False,
) -> Seq:
    consensus = ""
    con_len = alignment.get_alignment_length()

    for n in range(con_len):
        atom_dict: Counter[str] = Counter()
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


def _get_all_letters(alignment: Any) -> str:
    set_letters = set()
    for record in alignment:
        set_letters.update(record.seq)
    return "".join(sorted(set_letters))


def _pos_specific_score_matrix(
    alignment: Any,
    axis_seq: Seq | None = None,
    chars_to_ignore: list[str] | str | None = None,
) -> PSSM:
    all_letters = _get_all_letters(alignment)
    if not all_letters:
        raise ValueError("_get_all_letters returned empty string")

    if chars_to_ignore is None:
        chars_to_ignore = []
    elif isinstance(chars_to_ignore, str):
        chars_to_ignore = [chars_to_ignore]
    elif not isinstance(chars_to_ignore, list):
        raise TypeError("chars_to_ignore should be a list.")

    gap_char = "-"
    chars_to_ignore.append(gap_char)

    for char in chars_to_ignore:
        all_letters = all_letters.replace(char, "")

    if axis_seq:
        left_seq = axis_seq
        if len(axis_seq) != alignment.get_alignment_length():
            raise ValueError(
                "Axis sequence length does not equal the get_alignment_length"
            )
    else:
        left_seq = _dumb_consensus(alignment)

    pssm_info: list[tuple[str, dict[str, float]]] = []
    for residue_num in range(len(left_seq)):
        score_dict: dict[str, float] = dict.fromkeys(all_letters, 0)
        for record in alignment:
            try:
                this_residue = record.seq[residue_num]
            except IndexError:
                this_residue = None

            if this_residue and this_residue not in chars_to_ignore:
                weight = record.annotations.get("weight", 1.0)
                try:
                    score_dict[this_residue] += weight
                except KeyError:
                    raise ValueError(
                        "Residue %s not found" % this_residue
                    ) from None

        pssm_info.append((left_seq[residue_num], score_dict))

    return PSSM(pssm_info)


class PositionSpecificScoreMatrix(Alignment):
    def __init__(self, args: Any) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        if self.fasta is None:
            raise ValueError("fasta cannot be None")
        alignment = read_alignment_alignio(self.fasta)
        for seqrecord in alignment:
            seqrecord.seq._data = seqrecord.seq._data.upper().replace(b"?", b"-")
        consensus = _dumb_consensus(alignment)
        my_pssm = _pos_specific_score_matrix(
            alignment,
            consensus,
            chars_to_ignore=self.ambiguous_character,
        )
        print(vars(my_pssm))

    def process_args(self, args: Any) -> dict[str, Any]:
        return dict(
            fasta=args.fasta,
            ambiguous_character=args.ambiguous_character,
        )
