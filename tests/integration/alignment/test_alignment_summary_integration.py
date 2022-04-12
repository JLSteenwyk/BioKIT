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
    def test_alignment_summary_invalid_input(self, mocked_print):  # noqa
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Biokit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_alignment_summary_simple(self, mocked_print):
        expected_result = textwrap.dedent(
            """
            General Characteristics
            =======================
            5\tNumber of taxa
            6\tAlignment length
            3\tParsimony informative sites
            3\tVariable sites
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

    @patch("builtins.print")
    def test_alignment_summary_alias(self, mocked_print):
        expected_result = textwrap.dedent(
            """
            General Characteristics
            =======================
            5\tNumber of taxa
            6\tAlignment length
            3\tParsimony informative sites
            3\tVariable sites
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
            "aln_summary",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_alignment_summary_decode_test(self, mocked_print):
        expected_result = textwrap.dedent(
            """
            General Characteristics
            =======================
            8\tNumber of taxa
            17143\tAlignment length
            3767\tParsimony informative sites
            3767\tVariable sites
            13210\tConstant sites

            Character Frequencies
            =====================
            Y\t3745
            W\t1367
            V\t7925
            T\t6936
            S\t9949
            R\t7584
            Q\t5166
            P\t6397
            N\t4784
            M\t3104
            L\t13540
            K\t8252
            I\t6531
            H\t2777
            G\t7189
            F\t4945
            E\t9243
            D\t8140
            C\t1530
            A\t10442
            -\t7598"""
        )

        testargs = [
            "biokit",
            "aln_summary",
            f"{here.parent.parent.parent}/sample_files/busco10.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]
