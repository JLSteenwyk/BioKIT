import pytest

from mock import patch, call
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestNumberOfScaffolds(object):
    @patch("builtins.print")
    def test_number_of_scaffolds_invalid_input(self, mocked_print):  # noqa
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Biokit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_number_of_scaffolds_simple(self, mocked_print):
        expected_result = 17

        testargs = [
            "biokit",
            "number_of_scaffolds",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_genomic.fna",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_number_of_scaffolds_simple_alias0(self, mocked_print):
        expected_result = 17

        testargs = [
            "biokit",
            "number_of_large_scaffolds",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_genomic.fna",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_number_of_scaffolds_simple_alias1(self, mocked_print):
        expected_result = 17

        testargs = [
            "biokit",
            "num_of_scaffolds",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_genomic.fna",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_number_of_scaffolds_simple_alias2(self, mocked_print):
        expected_result = 17

        testargs = [
            "biokit",
            "number_of_contigs",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_genomic.fna",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_number_of_scaffolds_simple_alias3(self, mocked_print):
        expected_result = 17

        testargs = [
            "biokit",
            "num_of_cont",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_genomic.fna",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]
