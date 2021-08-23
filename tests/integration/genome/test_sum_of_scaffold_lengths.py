import pytest

from mock import patch, call
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestSumOfScaffoldLengths(object):
    @patch("builtins.print")
    def test_sum_of_scaffold_lengths_invalid_input(self, mocked_print):  # noqa
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Biokit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_sum_of_scaffold_lengths_simple(self, mocked_print):
        expected_result = 12157105

        testargs = [
            "biokit",
            "sum_of_scaffold_lengths",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_genomic.fna",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]
