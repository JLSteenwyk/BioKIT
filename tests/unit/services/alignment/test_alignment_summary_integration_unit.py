import pytest
from argparse import Namespace

from biokit.services.alignment.alignment_summary import AlignmentSummary


@pytest.fixture
def args():
    kwargs = dict(fasta="/some/path/to/file.fa")
    return Namespace(**kwargs)


class TestAlignmentSummary(object):
    def test_init_sets_fasta_file_path(self, args):
        aln = AlignmentSummary(args)
        assert aln.fasta == args.fasta
