import pytest

from mock import patch, call
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestFaidx(object):
    @patch("builtins.print")
    def test_faidx_invalid_input(self, mocked_print):  # noqa
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Biokit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_faidx(self, mocked_print):
        expected_result = ">1\nA?GTAT"
        testargs = [
            "biokit",
            "faidx",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            '-e',
            '1'
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_faidx_alias0(self, mocked_print):
        expected_result = ">1\nA?GTAT"
        testargs = [
            "biokit",
            "get_entry",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            '-e',
            '1'
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_faidx_alias1(self, mocked_print):
        expected_result = ">1\nA?GTAT"
        testargs = [
            "biokit",
            "ge",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            '-e',
            '1'
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]
