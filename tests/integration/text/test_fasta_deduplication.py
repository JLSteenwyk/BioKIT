import pytest

from mock import patch, call  # noqa
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestFastaDeduplication(object):
    @patch("builtins.print")
    def test_fasta_deduplication(self, mocked_print):
        testargs = [
            "biokit",
            "fasta_deduplication",
            f"{here.parent.parent.parent}/sample_files/duplicates.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/duplicates.dedup.fa", "r"
        ) as expected_out:
            expected_out = expected_out.read()

        with open(
            f"{here.parent.parent.parent}/sample_files/duplicates.fa.dedup.fa",
            "r",
        ) as output_file:
            output_file = output_file.read()

        assert expected_out == output_file

    @patch("builtins.print")
    def test_fasta_deduplication_alias(self, mocked_print):
        testargs = [
            "biokit",
            "dedup",
            f"{here.parent.parent.parent}/sample_files/duplicates.fa",
            "-o",
            f"{here.parent.parent.parent}/sample_files/duplicates.custom_dedup.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/duplicates.dedup.fa", "r"
        ) as expected_out:
            expected_out = expected_out.read()

        with open(
            f"{here.parent.parent.parent}/sample_files/duplicates.custom_dedup.fa",
            "r",
        ) as output_file:
            output_file = output_file.read()

        assert expected_out == output_file

    @patch("builtins.print")
    def test_fasta_deduplication_long_output_arg(self, mocked_print):
        testargs = [
            "biokit",
            "fasta_deduplication",
            f"{here.parent.parent.parent}/sample_files/duplicates.fa",
            "--output",
            f"{here.parent.parent.parent}/sample_files/duplicates.long_arg_dedup.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/duplicates.dedup.fa", "r"
        ) as expected_out:
            expected_out = expected_out.read()

        with open(
            f"{here.parent.parent.parent}/sample_files/duplicates.long_arg_dedup.fa",
            "r",
        ) as output_file:
            output_file = output_file.read()

        assert expected_out == output_file
