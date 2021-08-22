import pytest

from biokit.helpers.files import read_alignment_alignio


class TestFiles(object):
    def test_bad_alignment_file_path(self, mocker):
        in_file = ""
        mocker.patch("biokit.helpers.files.AlignIO.read", side_effect=ValueError())

        with pytest.raises(Exception) as excinfo:
            read_alignment_alignio(in_file)
        assert "Input file could not be read" in str(excinfo.value)
