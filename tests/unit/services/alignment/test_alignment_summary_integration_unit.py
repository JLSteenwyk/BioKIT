import pytest
from argparse import Namespace
from math import isclose
from mock import patch, call

from biokit.services.alignment.alignment_summary import AlignmentSummary
from biokit.services.alignment.base import Alignment


@pytest.fixture
def args():
    kwargs = dict(fasta="/some/path/to/file.fa")
    return Namespace(**kwargs)


class TestAlignmentSummary(object):
    def test_init_sets_fasta_file_path(self, args):
        aln = AlignmentSummary(args)
        assert aln.fasta == args.fasta
