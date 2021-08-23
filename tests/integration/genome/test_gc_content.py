import pytest

from mock import patch, call
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestGCContent(object):
    @patch("builtins.print")
    def test_gc_content_invalid_input(self, mocked_print):  # noqa
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Biokit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_gc_content_simple(self, mocked_print):
        expected_result = 0.2273

        testargs = [
            "biokit",
            "gc_content",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gc_content_verbose(self, mocked_print):
        expected_result_0 = "lcl|NC_001133.9_cds_NP_009332.1_1\t0.4959"
        expected_result_1 = "lcl|NC_001133.9_cds_NP_878038.1_2\t0.4123"
        expected_result_2 = "lcl|NC_001133.9_cds_NP_009333.1_3\t0.3608"

        testargs = [
            "biokit",
            "gc_content",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-v",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
        ]

    @pytest.mark.slow
    @patch("builtins.print")
    def test_gc_content_verbose_slow(self, mocked_print):
        expected_result_00 = "NC_001133.9\t0.3927"
        expected_result_01 = "NC_001134.8\t0.3834"
        expected_result_02 = "NC_001135.5\t0.3853"
        expected_result_03 = "NC_001136.10\t0.3791"
        expected_result_04 = "NC_001137.3\t0.3851"
        expected_result_05 = "NC_001138.5\t0.3873"
        expected_result_06 = "NC_001139.9\t0.3806"
        expected_result_07 = "NC_001140.6\t0.385"
        expected_result_08 = "NC_001141.2\t0.389"
        expected_result_09 = "NC_001142.9\t0.3837"
        expected_result_10 = "NC_001143.9\t0.3807"
        expected_result_11 = "NC_001144.5\t0.3848"
        expected_result_12 = "NC_001145.3\t0.382"
        expected_result_13 = "NC_001146.8\t0.3864"
        expected_result_14 = "NC_001147.6\t0.3816"
        expected_result_15 = "NC_001148.4\t0.3806"
        expected_result_16 = "NC_001224.1\t0.1711"

        testargs = [
            "biokit",
            "gc_content",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_genomic.fna",
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
        ]
