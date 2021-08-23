import pytest

from mock import patch, call
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestSequenceComplement(object):
    @patch("builtins.print")
    def test_character_frequency_invalid_input(self, mocked_print):  # noqa
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Biokit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_sequence_complement_simple(self, mocked_print):
        expected_result_0 = ">lcl|NC_001133.9_cds_NP_009332.1_1\nGACCAGTTTAATTGAAGTTAGCGGCGACCACAGCGACGGTAGCGACGATGACGAAGACGTTGGTGGTGAGATCGAGTTAGACTGCTTTCTCAGTTGAACCACCTTAACCCACAGATGCAGAGACTATAGTCTCGAGTGAATCGGGTTATGATGTACAAGGTTCGGCGGGTGGGTTGACTTTGGATGGGTCAGCTTCAACGACTTCGGCAAAAGTTGATGCCACTGAAGTGGTGGTACAACTGGCCATAACGAGGTCTGGTTCACTGGTCTTACTAGTGGCCACAAGGTACCATGAGGTCGTCTAATTTCGGTCGGTAGAGGTCACGAGATAGGTTCCTGCCATAGATGTGATAGCGTTTGATC"

        testargs = [
            "biokit",
            "sequence_complement",
            f"{here.parent.parent.parent}/sample_files/YAL068C_cds.fna",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result_0)]

    @patch("builtins.print")
    def test_sequence_complement_multiple_entries(self, mocked_print):
        expected_result_0 = ">200_S38\ntaccgactgtaggagtgcgtcgaggtctgaacggacctagtcgaacgttgtgttaagatgcgttgtgaaccaatagagtgttgtatggtgctgttacgggggtgttgtggtggtggtttacaggggctgcgtcgtggtcgggatcgtttctagtggttcttgagtagtagtggcggtcagggtcgtcggtagcgtttatttcaccccccacgtcgacgacaacgcccgttacgtagtgggggtgtccgcggaggagttgttgggcctcgacgcggtccctctcgtcatcttccacttctagggttggaaggagggcgcggtctgagcgggtcgtgcaaacgttcggccgtcgccctcgaacgcgcgctagagtaatagtttcttgtcgtctagctcatggaatagaggcacgaagggccctaaccgcggagactccgacttgttctttggtcttaggtcctggacctctggctcgaatctctgcagctcttcctcgcgcgacgctttcacgccctcaactttttcaactcctgagccaacctcctacaagaaccgcgacagcgacacccataggtgcccctaccaatgagagttttgact"
        expected_result_1 = ">203_S40\ntaccgactgtaggagtgcgtcgaggtctgaacggacctagtcgaacgttgtgttaagatgcgttgtgaaccaatagagtgttgtatggtgctgttacgggggtgttgtggtggtggtttacaggggctgcgtcgtggtcgggatcgtttctagtggttcttgagtagtagtggtggtcagggtcgtcggtagcgtttatttcaccccccacgtcgacgacaacgcccgttacgtagtgggggtgtccgcggaggagttgttgggcctcgacgcggtccctctcgtcatcttccacttctagggttggaaggagggcgcggtctgagcgggtcgtgcaaacgttcggccgtcgccctcgaacgcgcgctagagtaatagtttcttgtcgtctagctcatggaatagaggcacgaagggccctaaccgcggagactccgacttgttctttggtcttaggtcctggacctctggctcgaatctctgcagctcttcctcgcgcgacgctttcacgccctcaactttttcaactcctgagccaacctcctacaagaaccgcgacagcgacacccataggtgcccctaccaatgagagttttgact"

        testargs = [
            "biokit",
            "sequence_complement",
            f"{here.parent.parent.parent}/sample_files/EOG091N44MS.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
        ]

    @patch("builtins.print")
    def test_sequence_complement_simple_reverse(self, mocked_print):
        expected_result_0 = ">lcl|NC_001133.9_cds_NP_009332.1_1\nCTAGTTTGCGATAGTGTAGATACCGTCCTTGGATAGAGCACTGGAGATGGCTGGCTTTAATCTGCTGGAGTACCATGGAACACCGGTGATCATTCTGGTCACTTGGTCTGGAGCAATACCGGTCAACATGGTGGTGAAGTCACCGTAGTTGAAAACGGCTTCAGCAACTTCGACTGGGTAGGTTTCAGTTGGGTGGGCGGCTTGGAACATGTAGTATTGGGCTAAGTGAGCTCTGATATCAGAGACGTAGACACCCAATTCCACCAAGTTGACTCTTTCGTCAGATTGAGCTAGAGTGGTGGTTGCAGAAGCAGTAGCAGCGATGGCAGCGACACCAGCGGCGATTGAAGTTAATTTGACCAG"
        testargs = [
            "biokit",
            "sequence_complement",
            f"{here.parent.parent.parent}/sample_files/YAL068C_cds.fna",
            "-r",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result_0)]

    @patch("builtins.print")
    def test_sequence_complement_multiple_entries_reverse(self, mocked_print):
        expected_result_0 = ">200_S38\ntcagttttgagagtaaccatccccgtggatacccacagcgacagcgccaagaacatcctccaaccgagtcctcaactttttcaactcccgcactttcgcagcgcgctccttctcgacgtctctaagctcggtctccaggtcctggattctggtttcttgttcagcctcagaggcgccaatcccgggaagcacggagataaggtactcgatctgctgttctttgataatgagatcgcgcgcaagctcccgctgccggcttgcaaacgtgctgggcgagtctggcgcgggaggaaggttgggatcttcaccttctactgctctccctggcgcagctccgggttgttgaggaggcgcctgtgggggtgatgcattgcccgcaacagcagctgcaccccccactttatttgcgatggctgctgggactggcggtgatgatgagttcttggtgatctttgctagggctggtgctgcgtcggggacatttggtggtggtgttgtgggggcattgtcgtggtatgttgtgagataaccaagtgttgcgtagaattgtgttgcaagctgatccaggcaagtctggagctgcgtgaggatgtcagccat"
        expected_result_1 = ">203_S40\ntcagttttgagagtaaccatccccgtggatacccacagcgacagcgccaagaacatcctccaaccgagtcctcaactttttcaactcccgcactttcgcagcgcgctccttctcgacgtctctaagctcggtctccaggtcctggattctggtttcttgttcagcctcagaggcgccaatcccgggaagcacggagataaggtactcgatctgctgttctttgataatgagatcgcgcgcaagctcccgctgccggcttgcaaacgtgctgggcgagtctggcgcgggaggaaggttgggatcttcaccttctactgctctccctggcgcagctccgggttgttgaggaggcgcctgtgggggtgatgcattgcccgcaacagcagctgcaccccccactttatttgcgatggctgctgggactggtggtgatgatgagttcttggtgatctttgctagggctggtgctgcgtcggggacatttggtggtggtgttgtgggggcattgtcgtggtatgttgtgagataaccaagtgttgcgtagaattgtgttgcaagctgatccaggcaagtctggagctgcgtgaggatgtcagccat"
        testargs = [
            "biokit",
            "sequence_complement",
            f"{here.parent.parent.parent}/sample_files/EOG091N44MS.fa",
            "-r",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
        ]

    @patch("builtins.print")
    def test_sequence_complement_multiple_entries_alias(self, mocked_print):
        expected_result_0 = ">200_S38\ntaccgactgtaggagtgcgtcgaggtctgaacggacctagtcgaacgttgtgttaagatgcgttgtgaaccaatagagtgttgtatggtgctgttacgggggtgttgtggtggtggtttacaggggctgcgtcgtggtcgggatcgtttctagtggttcttgagtagtagtggcggtcagggtcgtcggtagcgtttatttcaccccccacgtcgacgacaacgcccgttacgtagtgggggtgtccgcggaggagttgttgggcctcgacgcggtccctctcgtcatcttccacttctagggttggaaggagggcgcggtctgagcgggtcgtgcaaacgttcggccgtcgccctcgaacgcgcgctagagtaatagtttcttgtcgtctagctcatggaatagaggcacgaagggccctaaccgcggagactccgacttgttctttggtcttaggtcctggacctctggctcgaatctctgcagctcttcctcgcgcgacgctttcacgccctcaactttttcaactcctgagccaacctcctacaagaaccgcgacagcgacacccataggtgcccctaccaatgagagttttgact"
        expected_result_1 = ">203_S40\ntaccgactgtaggagtgcgtcgaggtctgaacggacctagtcgaacgttgtgttaagatgcgttgtgaaccaatagagtgttgtatggtgctgttacgggggtgttgtggtggtggtttacaggggctgcgtcgtggtcgggatcgtttctagtggttcttgagtagtagtggtggtcagggtcgtcggtagcgtttatttcaccccccacgtcgacgacaacgcccgttacgtagtgggggtgtccgcggaggagttgttgggcctcgacgcggtccctctcgtcatcttccacttctagggttggaaggagggcgcggtctgagcgggtcgtgcaaacgttcggccgtcgccctcgaacgcgcgctagagtaatagtttcttgtcgtctagctcatggaatagaggcacgaagggccctaaccgcggagactccgacttgttctttggtcttaggtcctggacctctggctcgaatctctgcagctcttcctcgcgcgacgctttcacgccctcaactttttcaactcctgagccaacctcctacaagaaccgcgacagcgacacccataggtgcccctaccaatgagagttttgact"

        testargs = [
            "biokit",
            "seq_comp",
            f"{here.parent.parent.parent}/sample_files/EOG091N44MS.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
        ]
