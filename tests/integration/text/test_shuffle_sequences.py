import pytest

from mock import patch
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestShuffleSequences(object):
    @patch("builtins.print")
    def test_shuffle_sequences(self, mocked_print):
        testargs = [
            "biokit",
            "shuffle_sequences",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-o",
            f"{here.parent.parent.parent}/sample_files/simple.shuffled.fa",
            "-s",
            "42",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        out = Path(f"{here.parent.parent.parent}/sample_files/simple.shuffled.fa")
        assert out.exists()
        out.unlink()

    @patch("builtins.print")
    def test_shuffle_sequences_alias(self, mocked_print):
        testargs = [
            "biokit",
            "shuffle_seqs",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-o",
            f"{here.parent.parent.parent}/sample_files/simple.alias_shuffled.fa",
            "-s",
            "42",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        out = Path(f"{here.parent.parent.parent}/sample_files/simple.alias_shuffled.fa")
        assert out.exists()
        out.unlink()
