from argparse import Namespace

import pytest

from biokit.services.text.faidx import Faidx
from biokit.services.text.file_format_converter import FileFormatConverter
from biokit.services.text.rename_fasta_entries import RenameFastaEntries
from biokit.services.text.reorder_by_sequence_length import ReorderBySequenceLength


def test_file_format_converter_rejects_unknown_input_format(tmp_path, capsys):
    in_fasta = tmp_path / "in.fa"
    out_file = tmp_path / "out.fa"
    in_fasta.write_text(">a\nACGT\n")
    args = Namespace(
        input_file=str(in_fasta),
        input_file_format="unknown",
        output_file=str(out_file),
        output_file_format="fasta",
    )

    FileFormatConverter(args).run()
    assert "File format not acceptable" in capsys.readouterr().out
    assert not out_file.exists()


def test_reorder_by_sequence_length_handles_duplicate_ids(tmp_path):
    in_fasta = tmp_path / "in.fa"
    out_fasta = tmp_path / "out.fa"
    in_fasta.write_text(">dup\nAA\n>x\nAAAA\n>dup\nAAA\n")
    args = Namespace(fasta=str(in_fasta), output=str(out_fasta))

    ReorderBySequenceLength(args).run()
    out = out_fasta.read_text()
    assert out.count(">dup") == 2
    assert out.splitlines()[0] == ">x"


def test_rename_fasta_entries_ignores_comments_and_blank_lines(tmp_path):
    in_fasta = tmp_path / "in.fa"
    idmap = tmp_path / "idmap.txt"
    out_fasta = tmp_path / "out.fa"
    in_fasta.write_text(">a\nAA\n>b\nTT\n")
    idmap.write_text("# comment\n\na alpha\nb beta\n")
    args = Namespace(fasta=str(in_fasta), idmap=str(idmap), output=str(out_fasta))

    RenameFastaEntries(args).run()
    out = out_fasta.read_text()
    assert ">alpha" in out
    assert ">beta" in out


def test_faidx_raises_clear_error_for_missing_entry(tmp_path):
    fasta = tmp_path / "in.fa"
    fasta.write_text(">a\nAA\n")
    args = Namespace(fasta=str(fasta), entry="missing")

    with pytest.raises(Exception, match="Entry 'missing' not found"):
        Faidx(args).run()
