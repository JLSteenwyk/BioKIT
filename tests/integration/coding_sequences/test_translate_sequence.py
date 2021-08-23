import pytest

from mock import patch, call # noqa
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestRelativeSynonymousCodonUsage(object):
    @patch("builtins.print")
    def test_translate_sequence_invalid_input(self, mocked_print):  # noqa
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Biokit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @pytest.mark.slow
    @patch("builtins.print")
    def test_translate_sequence(self, mocked_print):
        input_file = f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.fna"

        testargs = [
            "biokit",
            "translate_sequence",
            input_file,
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.fna_trans_seq_default", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{input_file}.translated.fa", "r") as out_fa:
            out_fa_content = out_fa.read()

        assert expected_fa_content == out_fa_content

    @patch("builtins.print")
    def test_translate_sequence_tt1(self, mocked_print):
        input_file = f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna"

        testargs = [
            "biokit",
            "translate_sequence",
            input_file,
            "-tt",
            "1",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.small.fna.tt1", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{input_file}.translated.fa", "r") as out_fa:
            out_fa_content = out_fa.read()

        assert expected_fa_content == out_fa_content

    @patch("builtins.print")
    def test_translate_sequence_tt2(self, mocked_print):
        input_file = f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna"

        testargs = [
            "biokit",
            "translate_sequence",
            input_file,
            "-tt",
            "2",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.small.fna.tt2", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{input_file}.translated.fa", "r") as out_fa:
            out_fa_content = out_fa.read()

        assert expected_fa_content == out_fa_content

    @patch("builtins.print")
    def test_translate_sequence_tt3(self, mocked_print):
        input_file = f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna"

        testargs = [
            "biokit",
            "translate_sequence",
            input_file,
            "-tt",
            "3",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.small.fna.tt3", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{input_file}.translated.fa", "r") as out_fa:
            out_fa_content = out_fa.read()

        assert expected_fa_content == out_fa_content

    @patch("builtins.print")
    def test_translate_sequence_tt4(self, mocked_print):
        input_file = f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna"

        testargs = [
            "biokit",
            "translate_sequence",
            input_file,
            "-tt",
            "4",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.small.fna.tt4", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{input_file}.translated.fa", "r") as out_fa:
            out_fa_content = out_fa.read()

        assert expected_fa_content == out_fa_content

    @patch("builtins.print")
    def test_translate_sequence_tt5(self, mocked_print):
        input_file = f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna"

        testargs = [
            "biokit",
            "translate_sequence",
            input_file,
            "-tt",
            "5",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.small.fna.tt5", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{input_file}.translated.fa", "r") as out_fa:
            out_fa_content = out_fa.read()

        assert expected_fa_content == out_fa_content

    @patch("builtins.print")
    def test_translate_sequence_tt6(self, mocked_print):
        input_file = f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna"

        testargs = [
            "biokit",
            "translate_sequence",
            input_file,
            "-tt",
            "6",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.small.fna.tt6", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{input_file}.translated.fa", "r") as out_fa:
            out_fa_content = out_fa.read()

        assert expected_fa_content == out_fa_content

    @patch("builtins.print")
    def test_translate_sequence_tt9(self, mocked_print):
        input_file = f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna"

        testargs = [
            "biokit",
            "translate_sequence",
            input_file,
            "-tt",
            "9",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.small.fna.tt9", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{input_file}.translated.fa", "r") as out_fa:
            out_fa_content = out_fa.read()

        assert expected_fa_content == out_fa_content

    @patch("builtins.print")
    def test_translate_sequence_tt10(self, mocked_print):
        input_file = f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna"

        testargs = [
            "biokit",
            "translate_sequence",
            input_file,
            "-tt",
            "10",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.small.fna.tt10", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{input_file}.translated.fa", "r") as out_fa:
            out_fa_content = out_fa.read()

        assert expected_fa_content == out_fa_content

    @patch("builtins.print")
    def test_translate_sequence_tt11(self, mocked_print):
        input_file = f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna"

        testargs = [
            "biokit",
            "translate_sequence",
            input_file,
            "-tt",
            "11",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.small.fna.tt11", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{input_file}.translated.fa", "r") as out_fa:
            out_fa_content = out_fa.read()

        assert expected_fa_content == out_fa_content

    @patch("builtins.print")
    def test_translate_sequence_tt12(self, mocked_print):
        input_file = f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna"

        testargs = [
            "biokit",
            "translate_sequence",
            input_file,
            "-tt",
            "12",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.small.fna.tt12", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{input_file}.translated.fa", "r") as out_fa:
            out_fa_content = out_fa.read()

        assert expected_fa_content == out_fa_content

    @patch("builtins.print")
    def test_translate_sequence_tt13(self, mocked_print):
        input_file = f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna"

        testargs = [
            "biokit",
            "translate_sequence",
            input_file,
            "-tt",
            "13",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.small.fna.tt13", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{input_file}.translated.fa", "r") as out_fa:
            out_fa_content = out_fa.read()

        assert expected_fa_content == out_fa_content

    @patch("builtins.print")
    def test_translate_sequence_tt14(self, mocked_print):
        input_file = f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna"

        testargs = [
            "biokit",
            "translate_sequence",
            input_file,
            "-tt",
            "14",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.small.fna.tt14", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{input_file}.translated.fa", "r") as out_fa:
            out_fa_content = out_fa.read()

        assert expected_fa_content == out_fa_content

    @patch("builtins.print")
    def test_translate_sequence_tt16(self, mocked_print):
        input_file = f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna"

        testargs = [
            "biokit",
            "translate_sequence",
            input_file,
            "-tt",
            "16",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.small.fna.tt16", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{input_file}.translated.fa", "r") as out_fa:
            out_fa_content = out_fa.read()

        assert expected_fa_content == out_fa_content

    @patch("builtins.print")
    def test_translate_sequence_tt21(self, mocked_print):
        input_file = f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna"

        testargs = [
            "biokit",
            "translate_sequence",
            input_file,
            "-tt",
            "21",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.small.fna.tt21", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{input_file}.translated.fa", "r") as out_fa:
            out_fa_content = out_fa.read()

        assert expected_fa_content == out_fa_content

    @patch("builtins.print")
    def test_translate_sequence_tt22(self, mocked_print):
        input_file = f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna"

        testargs = [
            "biokit",
            "translate_sequence",
            input_file,
            "-tt",
            "22",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.small.fna.tt22", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{input_file}.translated.fa", "r") as out_fa:
            out_fa_content = out_fa.read()

        assert expected_fa_content == out_fa_content

    @patch("builtins.print")
    def test_translate_sequence_tt23(self, mocked_print):
        input_file = f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna"

        testargs = [
            "biokit",
            "translate_sequence",
            input_file,
            "-tt",
            "23",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.small.fna.tt23", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{input_file}.translated.fa", "r") as out_fa:
            out_fa_content = out_fa.read()

        assert expected_fa_content == out_fa_content

    @patch("builtins.print")
    def test_translate_sequence_tt24(self, mocked_print):
        input_file = f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna"

        testargs = [
            "biokit",
            "translate_sequence",
            input_file,
            "-tt",
            "24",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.small.fna.tt24", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{input_file}.translated.fa", "r") as out_fa:
            out_fa_content = out_fa.read()

        assert expected_fa_content == out_fa_content

    @patch("builtins.print")
    def test_translate_sequence_tt25(self, mocked_print):
        input_file = f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna"

        testargs = [
            "biokit",
            "translate_sequence",
            input_file,
            "-tt",
            "25",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.small.fna.tt25", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{input_file}.translated.fa", "r") as out_fa:
            out_fa_content = out_fa.read()

        assert expected_fa_content == out_fa_content

    @patch("builtins.print")
    def test_translate_sequence_tt26(self, mocked_print):
        input_file = f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna"

        testargs = [
            "biokit",
            "translate_sequence",
            input_file,
            "-tt",
            "26",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.small.fna.tt26", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{input_file}.translated.fa", "r") as out_fa:
            out_fa_content = out_fa.read()

        assert expected_fa_content == out_fa_content

    @patch("builtins.print")
    def test_translate_sequence_tt27(self, mocked_print):
        input_file = f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna"

        testargs = [
            "biokit",
            "translate_sequence",
            input_file,
            "-tt",
            "27",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.small.fna.tt27", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{input_file}.translated.fa", "r") as out_fa:
            out_fa_content = out_fa.read()

        assert expected_fa_content == out_fa_content

    @patch("builtins.print")
    def test_translate_sequence_tt28(self, mocked_print):
        input_file = f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna"

        testargs = [
            "biokit",
            "translate_sequence",
            input_file,
            "-tt",
            "28",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.small.fna.tt28", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{input_file}.translated.fa", "r") as out_fa:
            out_fa_content = out_fa.read()

        assert expected_fa_content == out_fa_content

    @patch("builtins.print")
    def test_translate_sequence_tt29(self, mocked_print):
        input_file = f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna"

        testargs = [
            "biokit",
            "translate_sequence",
            input_file,
            "-tt",
            "29",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.small.fna.tt29", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{input_file}.translated.fa", "r") as out_fa:
            out_fa_content = out_fa.read()

        assert expected_fa_content == out_fa_content

    @patch("builtins.print")
    def test_translate_sequence_tt30(self, mocked_print):
        input_file = f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna"

        testargs = [
            "biokit",
            "translate_sequence",
            input_file,
            "-tt",
            "30",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.small.fna.tt30", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{input_file}.translated.fa", "r") as out_fa:
            out_fa_content = out_fa.read()

        assert expected_fa_content == out_fa_content

    @patch("builtins.print")
    def test_translate_sequence_tt31(self, mocked_print):
        input_file = f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna"

        testargs = [
            "biokit",
            "translate_sequence",
            input_file,
            "-tt",
            "31",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.small.fna.tt31", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{input_file}.translated.fa", "r") as out_fa:
            out_fa_content = out_fa.read()

        assert expected_fa_content == out_fa_content

    @patch("builtins.print")
    def test_translate_sequence_tt33(self, mocked_print):
        input_file = f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna"

        testargs = [
            "biokit",
            "translate_sequence",
            input_file,
            "-tt",
            "33",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.small.fna.tt33", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{input_file}.translated.fa", "r") as out_fa:
            out_fa_content = out_fa.read()

        assert expected_fa_content == out_fa_content

    @patch("builtins.print")
    def test_translate_sequence_tt50(self, mocked_print):
        input_file = f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna"

        testargs = [
            "biokit",
            "translate_sequence",
            input_file,
            "-tt",
            "50",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.small.fna.tt50", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{input_file}.translated.fa", "r") as out_fa:
            out_fa_content = out_fa.read()

        assert expected_fa_content == out_fa_content

    @patch("builtins.print")
    def test_translate_sequence_tt_custom(self, mocked_print):
        input_file = f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna"

        testargs = [
            "biokit",
            "translate_sequence",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            '-tt',
            f"{here.parent.parent.parent}/sample_files/CUG_ala_code.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(f"{here.parent.parent}/expected/GCF_000146045.2_R64_cds_from_genomic.small.fna.tt50", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{input_file}.translated.fa", "r") as out_fa:
            out_fa_content = out_fa.read()

        assert expected_fa_content == out_fa_content
