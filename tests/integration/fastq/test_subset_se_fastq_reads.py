import pytest
import re

from mock import patch, call  # noqa
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestSubsetPEFastQReads(object):
    @patch("builtins.print")
    def test_subset_se_fastq_reads_invalid_input(self, mocked_print):  # noqa
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Biokit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_subset_se_fastq_reads(self, mocked_print):
        input_file_1 = (
            f"{here.parent.parent.parent}/sample_files/DRR284700_1_subset.fastq"
        )

        testargs = ["biokit", "subset_se_fastq_reads", input_file_1, "-s", "154"]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/DRR284700_1_subset_subset.fq", "r"
        ) as expected_fq_1:
            expected_fq_1 = expected_fq_1.read()

        output_file_1 = re.sub(".fastq$|.fq$", "_subset.fq", input_file_1)

        with open(output_file_1, "r") as output_fq_1:
            output_fq_1 = output_fq_1.read()

        assert expected_fq_1 == output_fq_1

    @patch("builtins.print")
    def test_subset_se_fastq_reads_alias(self, mocked_print):
        input_file_1 = (
            f"{here.parent.parent.parent}/sample_files/DRR284700_1_subset.fastq"
        )

        testargs = ["biokit", "subset_se_fastq", input_file_1, "-s", "154"]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/DRR284700_1_subset_subset.fq", "r"
        ) as expected_fq_1:
            expected_fq_1 = expected_fq_1.read()

        output_file_1 = re.sub(".fastq$|.fq$", "_subset.fq", input_file_1)

        with open(output_file_1, "r") as output_fq_1:
            output_fq_1 = output_fq_1.read()

        assert expected_fq_1 == output_fq_1
