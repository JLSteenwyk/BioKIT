import pytest

from mock import patch, call # noqa
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestRemoveShortSequences(object):
    @patch("builtins.print")
    def test_remove_short_sequences_invalid_input(self, mocked_print):  # noqa
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Biokit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_remove_short_sequences(self, mocked_print):
        testargs = [
            "biokit",
            "remove_short_sequences",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.fna",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.fna.long_seqs.fa", "r"
        ) as expected_out:
            expected_out = expected_out.read()

        with open(
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.fna.long_seqs.fa", "r"
        ) as output_file:
            output_file = output_file.read()

        assert expected_out == output_file

    @patch("builtins.print")
    def test_remove_short_sequences_long_arg(self, mocked_print):
        testargs = [
            "biokit",
            "remove_short_sequences",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.fna",
            "--output",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.fna.long_seqs.fa", "r"
        ) as expected_out:
            expected_out = expected_out.read()

        with open(
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa", "r"
        ) as output_file:
            output_file = output_file.read()

        assert expected_out == output_file

    @patch("builtins.print")
    def test_remove_short_sequences_custom_threshold(self, mocked_print):
        testargs = [
            "biokit",
            "remove_short_sequences",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.fna",
            "-t",
            "1500",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.1500_threshold.fna", "r"
        ) as expected_out:
            expected_out = expected_out.read()

        with open(
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.fna.long_seqs.fa", "r"
        ) as output_file:
            output_file = output_file.read()

        assert expected_out == output_file

    @patch("builtins.print")
    def test_remove_short_sequences_custom_threshold_long(self, mocked_print):
        testargs = [
            "biokit",
            "remove_short_sequences",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.fna",
            "--threshold",
            "1500",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.1500_threshold.fna", "r"
        ) as expected_out:
            expected_out = expected_out.read()

        with open(
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.fna.long_seqs.fa", "r"
        ) as output_file:
            output_file = output_file.read()

        assert expected_out == output_file
