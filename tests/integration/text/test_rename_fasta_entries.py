import pytest

from mock import patch, call # noqa
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestRenameFastaEntries(object):
    @patch("builtins.print")
    def test_rename_fasta_entries_invalid_input(self, mocked_print):  # noqa
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Biokit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_rename_fasta_entries(self, mocked_print):
        testargs = [
            "biokit",
            "rename_fasta_entries",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "--idmap",
            f"{here.parent.parent.parent}/sample_files/simple_fasta_idmap.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/simple.fa.renamed.fa", "r"
        ) as expected_out:
            expected_out = expected_out.read()

        with open(
            f"{here.parent.parent.parent}/sample_files/simple.fa.renamed.fa", "r"
        ) as output_file:
            output_file = output_file.read()

        assert expected_out == output_file

    @patch("builtins.print")
    def test_rename_fasta_entries_custom_out(self, mocked_print):
        testargs = [
            "biokit",
            "rename_fasta_entries",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "--idmap",
            f"{here.parent.parent.parent}/sample_files/simple_fasta_idmap.txt",
            "-o"
            f"{here.parent.parent.parent}/sample_files/simple.fa.custom_out",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/simple.fa.custom_out", "r"
        ) as expected_out:
            expected_out = expected_out.read()

        with open(
            f"{here.parent.parent.parent}/sample_files/simple.fa.custom_out", "r"
        ) as output_file:
            output_file = output_file.read()

        assert expected_out == output_file

    @patch("builtins.print")
    def test_rename_fasta_entries_alias(self, mocked_print):
        testargs = [
            "biokit",
            "rename_fasta",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "--idmap",
            f"{here.parent.parent.parent}/sample_files/simple_fasta_idmap.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/simple.fa.renamed.fa", "r"
        ) as expected_out:
            expected_out = expected_out.read()

        with open(
            f"{here.parent.parent.parent}/sample_files/simple.fa.renamed.fa", "r"
        ) as output_file:
            output_file = output_file.read()

        assert expected_out == output_file
