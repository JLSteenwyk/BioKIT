import pytest

from mock import patch, call # noqa
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestRemoveFastaEntry(object):
    @patch("builtins.print")
    def test_remove_fasta_entry_invalid_input(self, mocked_print):  # noqa
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Biokit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_remove_fasta_entry(self, mocked_print):
        testargs = [
            "biokit",
            "remove_fasta_entry",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-e",
            "1",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/simple_pruned.fa", "r"
        ) as expected_out:
            expected_out = expected_out.read()

        with open(
            f"{here.parent.parent.parent}/sample_files/simple.fa.pruned.fa", "r"
        ) as output_file:
            output_file = output_file.read()

        assert expected_out == output_file

    @patch("builtins.print")
    def test_remove_fasta_entry_long_arg(self, mocked_print):
        testargs = [
            "biokit",
            "remove_fasta_entry",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "--entry",
            "1",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/simple_pruned.fa", "r"
        ) as expected_out:
            expected_out = expected_out.read()

        with open(
            f"{here.parent.parent.parent}/sample_files/simple.fa.pruned.fa", "r"
        ) as output_file:
            output_file = output_file.read()

        assert expected_out == output_file

    @patch("builtins.print")
    def test_remove_fasta_entry_custom_out(self, mocked_print):
        testargs = [
            "biokit",
            "remove_fasta_entry",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "--entry",
            "1",
            "-o",
            f"{here.parent.parent.parent}/sample_files/simple.fa.custom_out.pruned.fa"
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/simple_pruned.fa", "r"
        ) as expected_out:
            expected_out = expected_out.read()

        with open(
            f"{here.parent.parent.parent}/sample_files/simple.fa.custom_out.pruned.fa", "r"
        ) as output_file:
            output_file = output_file.read()

        assert expected_out == output_file

    @patch("builtins.print")
    def test_remove_fasta_entry_custom_out_long_arg(self, mocked_print):
        testargs = [
            "biokit",
            "remove_fasta_entry",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "--entry",
            "1",
            "--output",
            f"{here.parent.parent.parent}/sample_files/simple.fa.custom_out.pruned.fa"
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/simple_pruned.fa", "r"
        ) as expected_out:
            expected_out = expected_out.read()

        with open(
            f"{here.parent.parent.parent}/sample_files/simple.fa.custom_out.pruned.fa", "r"
        ) as output_file:
            output_file = output_file.read()

        assert expected_out == output_file
