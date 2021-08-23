import pytest

from mock import patch, call # noqa
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestFileFormatConverter(object):
    @patch("builtins.print")
    def test_file_format_converter_invalid_input(self, mocked_print):  # noqa
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Biokit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_file_format_converter_fasta2clustal(self, mocked_print):
        testargs = [
            "biokit",
            "file_format_converter",
            "-i",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-iff",
            "fasta",
            "-o",
            f"{here.parent.parent.parent}/sample_files/simple.clustal",
            "-off",
            "clustal",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/simple.clustal", "r"
        ) as expected_out:
            expected_out = expected_out.read()

        with open(
            f"{here.parent.parent.parent}/sample_files/simple.clustal", "r"
        ) as output_file:
            output_file = output_file.read()

        assert expected_out == output_file

    @patch("builtins.print")
    def test_file_format_converter_fasta2maf(self, mocked_print):
        testargs = [
            "biokit",
            "file_format_converter",
            "-i",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-iff",
            "fasta",
            "-o",
            f"{here.parent.parent.parent}/sample_files/simple.maf",
            "-off",
            "maf",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/simple.maf", "r"
        ) as expected_out:
            expected_out = expected_out.read()

        with open(
            f"{here.parent.parent.parent}/sample_files/simple.maf", "r"
        ) as output_file:
            output_file = output_file.read()

        assert expected_out == output_file

    @patch("builtins.print")
    def test_file_format_converter_fasta2mauve(self, mocked_print):
        testargs = [
            "biokit",
            "file_format_converter",
            "-i",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-iff",
            "fasta",
            "-o",
            f"{here.parent.parent.parent}/sample_files/simple.mauve",
            "-off",
            "mauve",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/simple.mauve", "r"
        ) as expected_out:
            expected_out = expected_out.read()

        with open(
            f"{here.parent.parent.parent}/sample_files/simple.mauve", "r"
        ) as output_file:
            output_file = output_file.read()

        assert expected_out == output_file

    @patch("builtins.print")
    def test_file_format_converter_fasta2phylip(self, mocked_print):
        testargs = [
            "biokit",
            "file_format_converter",
            "-i",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-iff",
            "fasta",
            "-o",
            f"{here.parent.parent.parent}/sample_files/simple.phylip",
            "-off",
            "phylip",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/simple.phylip", "r"
        ) as expected_out:
            expected_out = expected_out.read()

        with open(
            f"{here.parent.parent.parent}/sample_files/simple.phylip", "r"
        ) as output_file:
            output_file = output_file.read()

        assert expected_out == output_file

    @patch("builtins.print")
    def test_file_format_converter_fasta2phylip_sequential(self, mocked_print):
        testargs = [
            "biokit",
            "file_format_converter",
            "-i",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-iff",
            "fasta",
            "-o",
            f"{here.parent.parent.parent}/sample_files/simple.phylip_sequential",
            "-off",
            "phylip_sequential",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/simple.phylip_sequential", "r"
        ) as expected_out:
            expected_out = expected_out.read()

        with open(
            f"{here.parent.parent.parent}/sample_files/simple.phylip_sequential", "r"
        ) as output_file:
            output_file = output_file.read()

        assert expected_out == output_file

    @patch("builtins.print")
    def test_file_format_converter_fasta2phylip_relaxed(self, mocked_print):
        testargs = [
            "biokit",
            "file_format_converter",
            "-i",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-iff",
            "fasta",
            "-o",
            f"{here.parent.parent.parent}/sample_files/simple.phylip_relaxed",
            "-off",
            "phylip_relaxed",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/simple.phylip_relaxed", "r"
        ) as expected_out:
            expected_out = expected_out.read()

        with open(
            f"{here.parent.parent.parent}/sample_files/simple.phylip_relaxed", "r"
        ) as output_file:
            output_file = output_file.read()

        assert expected_out == output_file

    @patch("builtins.print")
    def test_file_format_converter_fasta2stockholm(self, mocked_print):
        testargs = [
            "biokit",
            "file_format_converter",
            "-i",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-iff",
            "fasta",
            "-o",
            f"{here.parent.parent.parent}/sample_files/simple.stockholm",
            "-off",
            "stockholm",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/simple.stockholm", "r"
        ) as expected_out:
            expected_out = expected_out.read()

        with open(
            f"{here.parent.parent.parent}/sample_files/simple.stockholm", "r"
        ) as output_file:
            output_file = output_file.read()

        assert expected_out == output_file
