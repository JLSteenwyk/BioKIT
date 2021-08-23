import pytest
import re

from mock import patch, call  # noqa
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestTrimPEFastQ(object):
    @patch("builtins.print")
    def test_trim_pe_fastq_invalid_input(self, mocked_print):  # noqa
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Biokit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_trim_pe_fastq(self, mocked_print):
        input_file_1 = (
            f"{here.parent.parent.parent}/sample_files/DRR284700_1_subset.fastq"
        )
        input_file_2 = (
            f"{here.parent.parent.parent}/sample_files/DRR284700_2_subset.fastq"
        )

        testargs = [
            "biokit",
            "trim_pe_fastq",
            input_file_1,
            input_file_2,
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/DRR284700_1_subset_paired_trimmed.fq", "r"
        ) as expected_paired_fq_1, open(
            f"{here.parent.parent}/expected/DRR284700_1_subset_unpaired_trimmed.fq", "r"
        ) as expected_unpaired_fq_1, open(
            f"{here.parent.parent}/expected/DRR284700_2_subset_paired_trimmed.fq", "r"
        ) as expected_paired_fq_2, open(
            f"{here.parent.parent}/expected/DRR284700_2_subset_unpaired_trimmed.fq", "r"
        ) as expected_unpaired_fq_2:
            expected_paired_fq_1 = expected_paired_fq_1.read()
            expected_unpaired_fq_1 = expected_unpaired_fq_1.read()
            expected_paired_fq_2 = expected_paired_fq_2.read()
            expected_unpaired_fq_2 = expected_unpaired_fq_2.read()

        output_file_paired_1 = re.sub(
            ".fastq$|.fq$", "_paired_trimmed.fq", input_file_1
        )
        output_file_unpaired_1 = re.sub(
            ".fastq$|.fq$", "_unpaired_trimmed.fq", input_file_1
        )
        output_file_paired_2 = re.sub(
            ".fastq$|.fq$", "_paired_trimmed.fq", input_file_2
        )
        output_file_unpaired_2 = re.sub(
            ".fastq$|.fq$", "_unpaired_trimmed.fq", input_file_2
        )

        with open(output_file_paired_1, "r") as output_paired_fq_1, open(
            output_file_unpaired_1, "r"
        ) as output_unpaired_fq_1, open(
            output_file_paired_2, "r"
        ) as output_paired_fq_2, open(
            output_file_paired_2, "r"
        ) as output_paired_fq_2, open(
            output_file_unpaired_2, "r"
        ) as output_unpaired_fq_2:
            output_paired_fq_1 = output_paired_fq_1.read()
            output_unpaired_fq_1 = output_unpaired_fq_1.read()
            output_paired_fq_2 = output_paired_fq_2.read()
            output_unpaired_fq_2 = output_unpaired_fq_2.read()

        assert expected_paired_fq_1 == output_paired_fq_1
        assert expected_unpaired_fq_1 == output_unpaired_fq_1
        assert expected_paired_fq_2 == output_paired_fq_2
        assert expected_unpaired_fq_2 == output_unpaired_fq_2
