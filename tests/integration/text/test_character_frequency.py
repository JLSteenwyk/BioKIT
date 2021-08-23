import pytest

from mock import patch, call
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestCharacterFrequency(object):
    @patch("builtins.print")
    def test_character_frequency_invalid_input(self, mocked_print):  # noqa
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Biokit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @pytest.mark.slow
    @patch("builtins.print")
    def test_character_frequency_simple(self, mocked_print):
        expected_result_0 = "C\t0.1909"
        expected_result_1 = "A\t0.3098"
        expected_result_2 = "T\t0.3087"
        expected_result_3 = "G\t0.1906"

        testargs = [
            "biokit",
            "character_frequency",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_genomic.fna",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
            call(expected_result_3),
        ]

    @patch("builtins.print")
    def test_character_frequency_protein(self, mocked_print):
        expected_result_00 = "M\t0.0208"
        expected_result_01 = "R\t0.0444"
        expected_result_02 = "N\t0.0616"
        expected_result_03 = "E\t0.0653"
        expected_result_04 = "L\t0.0951"
        expected_result_05 = "Y\t0.0338"
        expected_result_06 = "Q\t0.0395"
        expected_result_07 = "W\t0.0104"
        expected_result_08 = "C\t0.0127"
        expected_result_09 = "V\t0.0556"
        expected_result_10 = "A\t0.0549"
        expected_result_11 = "S\t0.0898"
        expected_result_12 = "G\t0.0497"
        expected_result_13 = "K\t0.0734"
        expected_result_14 = "F\t0.0443"
        expected_result_15 = "T\t0.0591"
        expected_result_16 = "D\t0.0584"
        expected_result_17 = "I\t0.0656"
        expected_result_18 = "P\t0.0438"
        expected_result_19 = "H\t0.0217"

        testargs = [
            "biokit",
            "character_frequency",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_protein.faa",
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
