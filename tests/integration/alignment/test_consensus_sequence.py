import pytest

from mock import patch, call
from pathlib import Path
import sys

from biokit.biokit import Biokit


here = Path(__file__)


@pytest.mark.integration
class TestConsensusSequence(object):
    @patch("builtins.print")
    def test_consensus_sequence_invalid_input(self, mocked_print):  # noqa
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Biokit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_consensus_sequence_simple(self, mocked_print):
        expected_result = """>/Users/jlsteenwyk/Desktop/GITHUB/BioKIT/tests/sample_files/simple.fa.consensus\nANNTNN"""

        testargs = [
            "biokit",
            "consensus_sequence",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_consensus_sequence_simple_threshold(self, mocked_print):
        expected_result = """>/Users/jlsteenwyk/Desktop/GITHUB/BioKIT/tests/sample_files/simple.fa.consensus\nANGTTN"""

        testargs = [
            "biokit",
            "consensus_sequence",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            '-t',
            '.1',
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]
