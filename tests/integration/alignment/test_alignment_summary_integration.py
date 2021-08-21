import pytest

from mock import patch, call
from pathlib import Path
import sys
import textwrap

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestAlignmentSummary(object):
    @patch("builtins.print")
    def test_alignment_summary_invalid_input(self, mocked_print): # noqa
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Biokit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_alignment_summary(self, mocked_print):
        expected_result = textwrap.dedent(
            """
            General Characteristics
            =======================
            5\tNumber of taxa
            6\tAlignment length
            2\tParsimony informative sites
            2\tVariable sites
            1\tConstant sites

            Character Frequencies
            =====================
            T\t6
            G\t4
            C\t1
            A\t11
            -\t8"""
        )

        testargs = [
            "biokit",
            "alignment_summary",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]
