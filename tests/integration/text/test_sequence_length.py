import pytest

from mock import patch, call
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestSequenceLength(object):
    @patch("builtins.print")
    def test_character_frequency_invalid_input(self, mocked_print):  # noqa
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Biokit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_sequence_length_simple(self, mocked_print):
        expected_result_0 = "lcl|NC_001133.9_cds_NP_009332.1_1\t363"

        testargs = [
            "biokit",
            "sequence_length",
            f"{here.parent.parent.parent}/sample_files/YAL068C_cds.fna",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result_0)]

    @patch("builtins.print")
    def test_sequence_length_multiple_entries(self, mocked_print):
        expected_result_0 = "200_S38\t600"
        expected_result_1 = "203_S40\t600"

        testargs = [
            "biokit",
            "sequence_length",
            f"{here.parent.parent.parent}/sample_files/EOG091N44MS.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
        ]

    @patch("builtins.print")
    def test_sequence_length_simple_alias(self, mocked_print):
        expected_result_0 = "lcl|NC_001133.9_cds_NP_009332.1_1\t363"

        testargs = [
            "biokit",
            "seq_len",
            f"{here.parent.parent.parent}/sample_files/YAL068C_cds.fna",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result_0)]
