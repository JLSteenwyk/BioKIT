import pytest
import textwrap

from mock import patch, call
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestAlignmentRecoding(object):
    @patch("builtins.print")
    def test_alignment_recoding_invalid_input(self, mocked_print):  # noqa
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Biokit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_alignment_recoding_ry_nucleotide(self, mocked_print):
        expected_result_0 = textwrap.dedent(
            """>1\nR?RYRY"""
        )
        expected_result_1 = textwrap.dedent(
            """>2\nR-R-RY"""
        )
        expected_result_2 = textwrap.dedent(
            """>3\nR-R-YR"""
        )
        expected_result_3 = textwrap.dedent(
            """>4\nRRR-YR"""
        )
        expected_result_4 = textwrap.dedent(
            """>5\nRYR-Y-"""
        )

        testargs = [
            "biokit",
            "alignment_recoding",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-c",
            "RY-nucleotide",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
            call(expected_result_3),
            call(expected_result_4),
        ]

    @patch("builtins.print")
    def test_alignment_recoding_dayhoff_6(self, mocked_print):
        expected_result_0 = textwrap.dedent(
            """>200_S38\n301330131053113001440030430042110000000--130100003023021000003000301230000030010-000100---0---------------11000---00-02031011013000010000400212130213332111314330330030001011102311310132131212002321322320231133003030320104011---"""
        )
        expected_result_1 = textwrap.dedent(
            """>203_S40\n301330131053113001440030430042110000000--130100003023021000003000301230000030010-000100---0---------------11000---00-02031011013000010000400212130213332111314330330030001011102311310132131212002321322320231133003030320104011---"""
        )

        testargs = [
            "biokit",
            "alignment_recoding",
            f"{here.parent.parent.parent}/sample_files/EOG091N44MS.fa.mafft",
            "-c",
            "Dayhoff-6",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
        ]

    @patch("builtins.print")
    def test_alignment_recoding_sandr_6(self, mocked_print):
        expected_result_0 = textwrap.dedent(
            """>200_S38\n301330232043123002550031530055110000000--130100003023021000003000301231100030110-000200---0---------------22010---00-12031111013000010000500222130213332122315330330131001012102321310132131212002321322320231133103031351115021---"""
        )
        expected_result_1 = textwrap.dedent(
            """>203_S40\n301330232043123002550031530055110000000--130100003023021000003000301231100030110-000200---0---------------22010---00-12031111013000010000500222130213332122315330330131001012102321310132131212002321322320231133103031351115021---"""
        )

        testargs = [
            "biokit",
            "alignment_recoding",
            f"{here.parent.parent.parent}/sample_files/EOG091N44MS.fa.mafft",
            "-c",
            "SandR-6",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
        ]

    @patch("builtins.print")
    def test_alignment_recoding_KGB_6(self, mocked_print):
        expected_result_0 = textwrap.dedent(
            """>200_S38\n201221121152112011440120421141110011000--150100002012111000005000201150000050010-000100---0---------------11000---00-01051011012000010001400111120112221111214220520020001011111211211121151111001511211211121152005050210104011---"""
        )
        expected_result_1 = textwrap.dedent(
            """>203_S40\n201221121152112011440120421141110011000--150100002012111000005000201150000050010-000100---0---------------11000---00-01051011012000010001400111120112221111214220520020001011111211211121151111001511211211121152005050210104011---"""
        )

        testargs = [
            "biokit",
            "alignment_recoding",
            f"{here.parent.parent.parent}/sample_files/EOG091N44MS.fa.mafft",
            "-c",
            "KGB-6",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
        ]

    @patch("builtins.print")
    def test_alignment_recoding_alias0(self, mocked_print):
        expected_result_0 = textwrap.dedent(
            """>1\nR?RYRY"""
        )
        expected_result_1 = textwrap.dedent(
            """>2\nR-R-RY"""
        )
        expected_result_2 = textwrap.dedent(
            """>3\nR-R-YR"""
        )
        expected_result_3 = textwrap.dedent(
            """>4\nRRR-YR"""
        )
        expected_result_4 = textwrap.dedent(
            """>5\nRYR-Y-"""
        )

        testargs = [
            "biokit",
            "aln_recoding",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-c",
            "RY-nucleotide",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
            call(expected_result_3),
            call(expected_result_4),
        ]

    @patch("builtins.print")
    def test_alignment_recoding_alias1(self, mocked_print):
        expected_result_0 = textwrap.dedent(
            """>1\nR?RYRY"""
        )
        expected_result_1 = textwrap.dedent(
            """>2\nR-R-RY"""
        )
        expected_result_2 = textwrap.dedent(
            """>3\nR-R-YR"""
        )
        expected_result_3 = textwrap.dedent(
            """>4\nRRR-YR"""
        )
        expected_result_4 = textwrap.dedent(
            """>5\nRYR-Y-"""
        )

        testargs = [
            "biokit",
            "recode",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-c",
            "RY-nucleotide",

        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
            call(expected_result_3),
            call(expected_result_4),
        ]
    
    @patch("builtins.print")
    def test_alignment_recoding_no_recoding_table(self, mocked_print):
        expected_call = textwrap.dedent(
            """Please specify a recoding table"""
        )

        testargs = [
            "biokit",
            "recode",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]

        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Biokit()

        assert pytest_wrapped_e.type == SystemExit
        mocked_print.assert_has_calls([
            call(expected_call),
        ])
