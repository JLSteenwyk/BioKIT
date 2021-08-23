import pytest

from mock import patch, call
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestFastQReadLengths(object):
    @patch("builtins.print")
    def test_fastq_read_lengths_invalid_input(self, mocked_print):  # noqa
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Biokit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_fastq_read_read_lengths_simple(self, mocked_print):
        expected_result = "292.15 +/- 20.4019"

        testargs = [
            "biokit",
            "fastq_read_lengths",
            f"{here.parent.parent.parent}/sample_files/DRR284700_1_subset.fastq",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_fastq_read_read_lengths_simple_verbose(self, mocked_print):
        expected_result_00 = 299
        expected_result_01 = 234
        expected_result_02 = 299
        expected_result_03 = 298
        expected_result_04 = 301
        expected_result_05 = 261
        expected_result_06 = 300
        expected_result_07 = 301
        expected_result_08 = 301
        expected_result_09 = 300
        expected_result_10 = 301
        expected_result_11 = 301
        expected_result_12 = 300
        expected_result_13 = 243
        expected_result_14 = 301
        expected_result_15 = 300
        expected_result_16 = 301
        expected_result_17 = 301
        expected_result_18 = 300
        expected_result_19 = 301

        testargs = [
            "biokit",
            "fastq_read_lengths",
            f"{here.parent.parent.parent}/sample_files/DRR284700_1_subset.fastq",
            "-v",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [
            call(expected_result_00),
            call(expected_result_01),
            call(expected_result_02),
            call(expected_result_03),
            call(expected_result_04),
            call(expected_result_05),
            call(expected_result_06),
            call(expected_result_07),
            call(expected_result_08),
            call(expected_result_09),
            call(expected_result_10),
            call(expected_result_11),
            call(expected_result_12),
            call(expected_result_13),
            call(expected_result_14),
            call(expected_result_15),
            call(expected_result_16),
            call(expected_result_17),
            call(expected_result_18),
            call(expected_result_19),
        ]

    @patch("builtins.print")
    def test_fastq_read_read_lengths_simple_alias(self, mocked_print):
        expected_result = "292.15 +/- 20.4019"

        testargs = [
            "biokit",
            "fastq_read_lens",
            f"{here.parent.parent.parent}/sample_files/DRR284700_1_subset.fastq",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]
