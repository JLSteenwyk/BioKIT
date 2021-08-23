import pytest

from mock import patch, call
from pathlib import Path
import sys

from biokit.biokit import Biokit


here = Path(__file__)


@pytest.mark.integration
class TestConsensusSequence(object):
    @patch("builtins.print")
    def test_pssm_invalid_input(self, mocked_print):  # noqa
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Biokit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_pssm_simple(self, mocked_print):
        expected_result = {
            "pssm": [
                ("A", {"-": 0, "A": 5.0, "C": 0, "G": 0, "T": 0}),
                ("X", {"-": 3.0, "A": 0, "C": 1.0, "G": 1.0, "T": 0}),
                ("X", {"-": 0, "A": 2.0, "C": 0, "G": 3.0, "T": 0}),
                ("T", {"-": 4.0, "A": 0, "C": 0, "G": 0, "T": 1.0}),
                ("X", {"-": 0, "A": 2.0, "C": 0, "G": 0, "T": 3.0}),
                ("X", {"-": 1.0, "A": 2.0, "C": 0, "G": 0, "T": 2.0}),
            ]
        }

        testargs = [
            "biokit",
            "position_specific_score_matrix",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_pssm_simple_alias(self, mocked_print):
        expected_result = {
            "pssm": [
                ("A", {"-": 0, "A": 5.0, "C": 0, "G": 0, "T": 0}),
                ("X", {"-": 3.0, "A": 0, "C": 1.0, "G": 1.0, "T": 0}),
                ("X", {"-": 0, "A": 2.0, "C": 0, "G": 3.0, "T": 0}),
                ("T", {"-": 4.0, "A": 0, "C": 0, "G": 0, "T": 1.0}),
                ("X", {"-": 0, "A": 2.0, "C": 0, "G": 0, "T": 3.0}),
                ("X", {"-": 1.0, "A": 2.0, "C": 0, "G": 0, "T": 2.0}),
            ]
        }

        testargs = [
            "biokit",
            "pssm",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]
