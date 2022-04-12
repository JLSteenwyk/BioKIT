import pytest
import re

from mock import patch, call  # noqa
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestTrimPEAdaptersFastQ(object):
    @patch("builtins.print")
    def test_trim_pe_adapters_fastq_invalid_input(self, mocked_print):  # noqa
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Biokit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_trim_pe_adapters_fastq(self, mocked_print):
        input_file_1 = (
            f"{here.parent.parent.parent}/sample_files/DRR284700_1_subset.fastq"
        )
        input_file_2 = (
            f"{here.parent.parent.parent}/sample_files/DRR284700_2_subset.fastq"
        )

        testargs = [
            "biokit",
            "trim_pe_adapters_fastq",
            input_file_1,
            input_file_2,
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/DRR284700_1_subset_paired_adapter_trimmed.fq", "r"
        ) as expected_paired_fq_1, open(
            f"{here.parent.parent}/expected/DRR284700_2_subset_paired_adapter_trimmed.fq", "r"
        ) as expected_paired_fq_2:
            expected_paired_fq_1 = expected_paired_fq_1.read()
            expected_paired_fq_2 = expected_paired_fq_2.read()

        output_file_paired_1 = re.sub(
            ".fastq$|.fq$", "_paired_adapter_trimmed.fq", input_file_1
        )
        output_file_paired_2 = re.sub(
            ".fastq$|.fq$", "_paired_adapter_trimmed.fq", input_file_2
        )

        with open(output_file_paired_1, "r") as output_paired_fq_1, open(
            output_file_paired_2, "r"
        ) as output_paired_fq_2:
            output_paired_fq_1 = output_paired_fq_1.read()
            output_paired_fq_2 = output_paired_fq_2.read()

        assert expected_paired_fq_1 == output_paired_fq_1
        assert expected_paired_fq_2 == output_paired_fq_2

    @patch("builtins.print")
    def test_trim_pe_adapters_fastq_a0(self, mocked_print):
        input_file_1 = (
            f"{here.parent.parent.parent}/sample_files/DRR284700_1_subset.fastq"
        )
        input_file_2 = (
            f"{here.parent.parent.parent}/sample_files/DRR284700_2_subset.fastq"
        )

        testargs = [
            "biokit",
            "trim_pe_adapters_fastq",
            input_file_1,
            input_file_2,
            "-a",
            "NexteraPE-PE"
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/DRR284700_1_subset_paired_adapter_trimmed_NexteraPE-PE.fq", "r"
        ) as expected_paired_fq_1, open(
            f"{here.parent.parent}/expected/DRR284700_2_subset_paired_adapter_trimmed_NexteraPE-PE.fq", "r"
        ) as expected_paired_fq_2:
            expected_paired_fq_1 = expected_paired_fq_1.read()
            expected_paired_fq_2 = expected_paired_fq_2.read()

        output_file_paired_1 = re.sub(
            ".fastq$|.fq$", "_paired_adapter_trimmed.fq", input_file_1
        )
        output_file_paired_2 = re.sub(
            ".fastq$|.fq$", "_paired_adapter_trimmed.fq", input_file_2
        )

        with open(output_file_paired_1, "r") as output_paired_fq_1, open(
            output_file_paired_2, "r"
        ) as output_paired_fq_2:
            output_paired_fq_1 = output_paired_fq_1.read()
            output_paired_fq_2 = output_paired_fq_2.read()

        assert expected_paired_fq_1 == output_paired_fq_1
        assert expected_paired_fq_2 == output_paired_fq_2

    @patch("builtins.print")
    def test_trim_pe_adapters_fastq_a1(self, mocked_print):
        input_file_1 = (
            f"{here.parent.parent.parent}/sample_files/DRR284700_1_subset.fastq"
        )
        input_file_2 = (
            f"{here.parent.parent.parent}/sample_files/DRR284700_2_subset.fastq"
        )

        testargs = [
            "biokit",
            "trim_pe_adapters_fastq",
            input_file_1,
            input_file_2,
            "-a",
            "TruSeq2-PE"
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/DRR284700_1_subset_paired_adapter_trimmed_TruSeq2-PE.fq", "r"
        ) as expected_paired_fq_1, open(
            f"{here.parent.parent}/expected/DRR284700_2_subset_paired_adapter_trimmed_TruSeq2-PE.fq", "r"
        ) as expected_paired_fq_2:
            expected_paired_fq_1 = expected_paired_fq_1.read()
            expected_paired_fq_2 = expected_paired_fq_2.read()

        output_file_paired_1 = re.sub(
            ".fastq$|.fq$", "_paired_adapter_trimmed.fq", input_file_1
        )
        output_file_paired_2 = re.sub(
            ".fastq$|.fq$", "_paired_adapter_trimmed.fq", input_file_2
        )

        with open(output_file_paired_1, "r") as output_paired_fq_1, open(
            output_file_paired_2, "r"
        ) as output_paired_fq_2:
            output_paired_fq_1 = output_paired_fq_1.read()
            output_paired_fq_2 = output_paired_fq_2.read()

        assert expected_paired_fq_1 == output_paired_fq_1
        assert expected_paired_fq_2 == output_paired_fq_2

    @patch("builtins.print")
    def test_trim_pe_adapters_fastq_a2(self, mocked_print):
        input_file_1 = (
            f"{here.parent.parent.parent}/sample_files/DRR284700_1_subset.fastq"
        )
        input_file_2 = (
            f"{here.parent.parent.parent}/sample_files/DRR284700_2_subset.fastq"
        )

        testargs = [
            "biokit",
            "trim_pe_adapters_fastq",
            input_file_1,
            input_file_2,
            "-a",
            "TruSeq3-PE-2"
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/DRR284700_1_subset_paired_adapter_trimmed_TruSeq3-PE-2.fq", "r"
        ) as expected_paired_fq_1, open(
            f"{here.parent.parent}/expected/DRR284700_2_subset_paired_adapter_trimmed_TruSeq3-PE-2.fq", "r"
        ) as expected_paired_fq_2:
            expected_paired_fq_1 = expected_paired_fq_1.read()
            expected_paired_fq_2 = expected_paired_fq_2.read()

        output_file_paired_1 = re.sub(
            ".fastq$|.fq$", "_paired_adapter_trimmed.fq", input_file_1
        )
        output_file_paired_2 = re.sub(
            ".fastq$|.fq$", "_paired_adapter_trimmed.fq", input_file_2
        )

        with open(output_file_paired_1, "r") as output_paired_fq_1, open(
            output_file_paired_2, "r"
        ) as output_paired_fq_2:
            output_paired_fq_1 = output_paired_fq_1.read()
            output_paired_fq_2 = output_paired_fq_2.read()

        assert expected_paired_fq_1 == output_paired_fq_1
        assert expected_paired_fq_2 == output_paired_fq_2

    @patch("builtins.print")
    def test_trim_pe_adapters_fastq_a3(self, mocked_print):
        input_file_1 = (
            f"{here.parent.parent.parent}/sample_files/DRR284700_1_subset.fastq"
        )
        input_file_2 = (
            f"{here.parent.parent.parent}/sample_files/DRR284700_2_subset.fastq"
        )

        testargs = [
            "biokit",
            "trim_pe_adapters_fastq",
            input_file_1,
            input_file_2,
            "-a",
            "TruSeq3-PE"
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/DRR284700_1_subset_paired_adapter_trimmed_TruSeq3-PE.fq", "r"
        ) as expected_paired_fq_1, open(
            f"{here.parent.parent}/expected/DRR284700_2_subset_paired_adapter_trimmed_TruSeq3-PE.fq", "r"
        ) as expected_paired_fq_2:
            expected_paired_fq_1 = expected_paired_fq_1.read()
            expected_paired_fq_2 = expected_paired_fq_2.read()

        output_file_paired_1 = re.sub(
            ".fastq$|.fq$", "_paired_adapter_trimmed.fq", input_file_1
        )
        output_file_paired_2 = re.sub(
            ".fastq$|.fq$", "_paired_adapter_trimmed.fq", input_file_2
        )

        with open(output_file_paired_1, "r") as output_paired_fq_1, open(
            output_file_paired_2, "r"
        ) as output_paired_fq_2:
            output_paired_fq_1 = output_paired_fq_1.read()
            output_paired_fq_2 = output_paired_fq_2.read()

        assert expected_paired_fq_1 == output_paired_fq_1
        assert expected_paired_fq_2 == output_paired_fq_2

    @patch("builtins.print")
    def test_trim_pe_adapters_fastq_a4(self, mocked_print):
        input_file_1 = (
            f"{here.parent.parent.parent}/sample_files/SRR519926.1_1_subset_1.fastq"
        )
        input_file_2 = (
            f"{here.parent.parent.parent}/sample_files/SRR519926.1_1_subset_2.fastq"
        )

        testargs = [
            "biokit",
            "trim_pe_adapters_fastq",
            input_file_1,
            input_file_2,
            "-a",
            "TruSeq3-SE"
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/SRR519926.1_1_subset_1_paired_adapter_trimmed_TruSeq3-SE.fq", "r"
        ) as expected_paired_fq_1, open(
            f"{here.parent.parent}/expected/SRR519926.1_1_subset_2_paired_adapter_trimmed_TruSeq3-SE.fq", "r"
        ) as expected_paired_fq_2:
            expected_paired_fq_1 = expected_paired_fq_1.read()
            expected_paired_fq_2 = expected_paired_fq_2.read()

        output_file_paired_1 = re.sub(
            ".fastq$|.fq$", "_paired_adapter_trimmed.fq", input_file_1
        )
        output_file_paired_2 = re.sub(
            ".fastq$|.fq$", "_paired_adapter_trimmed.fq", input_file_2
        )

        with open(output_file_paired_1, "r") as output_paired_fq_1, open(
            output_file_paired_2, "r"
        ) as output_paired_fq_2:
            output_paired_fq_1 = output_paired_fq_1.read()
            output_paired_fq_2 = output_paired_fq_2.read()

        assert expected_paired_fq_1 == output_paired_fq_1
        assert expected_paired_fq_2 == output_paired_fq_2
