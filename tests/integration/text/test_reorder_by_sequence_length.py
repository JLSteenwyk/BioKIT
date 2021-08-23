import pytest

from mock import patch, call # noqa
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestReorderBySequenceLength(object):
    @patch("builtins.print")
    def test_reorder_by_sequence_length_invalid_input(self, mocked_print):  # noqa
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Biokit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_reorder_by_sequence_length(self, mocked_print):
        testargs = [
            "biokit",
            "reorder_by_sequence_length",
            f"{here.parent.parent.parent}/sample_files/EOG091N44MS.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/EOG091N44MS.fa.reordered.fa", "r"
        ) as expected_out:
            expected_out = expected_out.read()

        with open(
            f"{here.parent.parent.parent}/sample_files/EOG091N44MS.fa.reordered.fa", "r"
        ) as output_file:
            output_file = output_file.read()

        assert expected_out == output_file
