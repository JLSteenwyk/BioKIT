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
    def test_gc_content_third_position_invalid_input(self, mocked_print):  # noqa
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Biokit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @pytest.mark.slow
    @patch("builtins.print")
    def test_gc_content_third_position_simple(self, mocked_print):
        expected_result = "0.1222"

        testargs = [
            "biokit",
            "gc_content_third_position",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.fna",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gc_content_third_position_simple_verbose(self, mocked_print):
        expected_result_0 = textwrap.dedent(
            """lcl|NC_001133.9_cds_NP_009332.1_1\t0.4876"""
        )
        expected_result_1 = textwrap.dedent(
            """lcl|NC_001133.9_cds_NP_878038.1_2\t0.4211"""
        )
        expected_result_2 = textwrap.dedent(
            """lcl|NC_001133.9_cds_NP_009333.1_3\t0.3956"""
        )

        testargs = [
            "biokit",
            "gc_content_first_position",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            '-v',
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
        ]
