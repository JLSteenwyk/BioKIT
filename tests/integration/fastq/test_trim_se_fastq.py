import pytest
import re

from mock import patch, call  # noqa
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestTrimSEFastQ(object):
    @patch("builtins.print")
    def test_trim_se_fastq_invalid_input(self, mocked_print):  # noqa
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Biokit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_trim_se_fastq(self, mocked_print):
        input_file_1 = (
            f"{here.parent.parent.parent}/sample_files/DRR284700_1_subset.fastq"
        )

        testargs = [
            "biokit",
            "trim_se_fastq",
            input_file_1,
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        with open(
            f"{here.parent.parent}/expected/DRR284700_1_subset_trimmed.fq", "r"
        ) as expected_trimmed:
            expected_trimmed = expected_trimmed.read()

        output_file_trimmed = re.sub(".fastq$|.fq$", "_trimmed.fq", input_file_1)

        with open(output_file_trimmed, "r") as output_file_trimmed:
            output_file_trimmed = output_file_trimmed.read()

        assert expected_trimmed == output_file_trimmed
