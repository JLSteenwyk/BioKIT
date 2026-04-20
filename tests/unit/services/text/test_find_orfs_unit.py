import json
from argparse import Namespace

import pytest

from biokit.services.text.find_orfs import FindOrfs


def _args(fasta, min_length=100, translation_table=1, extract=False,
          protein=False, format=None):
    return Namespace(
        fasta=fasta,
        min_length=min_length,
        translation_table=translation_table,
        extract=extract,
        protein=protein,
        format=format,
    )


def _write_fasta(tmp_path, seqs):
    fa = tmp_path / "seqs.fa"
    with open(fa, "w") as fh:
        for name, seq in seqs:
            fh.write(f">{name}\n{seq}\n")
    return str(fa)


# A clean 100-aa forward ORF: ATG + 99 Ala codons (GCT) + TAA = 303 nt
CLEAN_ORF_FWD = "ATG" + "GCT" * 99 + "TAA"


def test_find_orfs_simple_forward(tmp_path, capsys):
    # Prepend 6 nt noise so ORF starts at position 7
    seq = "AAAAAA" + CLEAN_ORF_FWD + "CCCCCC"
    fa = _write_fasta(tmp_path, [("chr", seq)])
    FindOrfs(_args(fa, min_length=50)).run()
    lines = capsys.readouterr().out.strip().split("\n")
    assert lines[0] == "id\torf_id\tframe\tstart\tstop\tlength_nt\tlength_aa"
    assert len(lines) == 2
    fields = lines[1].split("\t")
    assert fields[0] == "chr"
    assert fields[2] == "+1"
    assert fields[3] == "7"
    assert fields[4] == "309"
    assert fields[5] == "303"
    assert fields[6] == "100"


def test_find_orfs_min_length_filter(tmp_path, capsys):
    # 50-aa ORF, filter requires 100
    short_orf = "ATG" + "GCT" * 49 + "TAA"
    fa = _write_fasta(tmp_path, [("chr", "AAA" + short_orf + "CCC")])
    FindOrfs(_args(fa, min_length=100)).run()
    out = capsys.readouterr().out.strip()
    # Only the header; no matching ORFs
    assert out == "id\torf_id\tframe\tstart\tstop\tlength_nt\tlength_aa"


def test_find_orfs_reports_reverse_strand(tmp_path, capsys):
    # Put the clean ORF on the reverse strand by reverse-complementing
    from Bio.Seq import Seq
    rc_of_orf = str(Seq(CLEAN_ORF_FWD).reverse_complement())
    seq = "AAAAAA" + rc_of_orf + "CCCCCC"
    fa = _write_fasta(tmp_path, [("chr", seq)])
    FindOrfs(_args(fa, min_length=50)).run()
    lines = capsys.readouterr().out.strip().split("\n")
    assert len(lines) == 2
    fields = lines[1].split("\t")
    # Frame should be a reverse-strand frame
    assert fields[2].startswith("-")
    assert fields[6] == "100"


def test_find_orfs_multiple_orfs_in_one_record(tmp_path, capsys):
    # Two forward ORFs separated by a short region
    seq = "AAA" + CLEAN_ORF_FWD + "TTT" + CLEAN_ORF_FWD + "CCC"
    fa = _write_fasta(tmp_path, [("chr", seq)])
    FindOrfs(_args(fa, min_length=50)).run()
    lines = capsys.readouterr().out.strip().split("\n")
    # Header + at least 2 ORFs
    assert len(lines) >= 3


def test_find_orfs_extract_nt(tmp_path, capsys):
    seq = "AAAAAA" + CLEAN_ORF_FWD + "CCCCCC"
    fa = _write_fasta(tmp_path, [("chr", seq)])
    FindOrfs(_args(fa, min_length=50, extract=True)).run()
    out = capsys.readouterr().out
    lines = out.strip().split("\n")
    # First line is FASTA header
    assert lines[0].startswith(">")
    assert "source=chr" in lines[0]
    assert "frame=+1" in lines[0]
    # Second line is the ORF nucleotide sequence
    assert lines[1].startswith("ATG")
    assert lines[1].endswith("TAA")
    assert len(lines[1]) == 303


def test_find_orfs_extract_protein(tmp_path, capsys):
    seq = "AAAAAA" + CLEAN_ORF_FWD + "CCCCCC"
    fa = _write_fasta(tmp_path, [("chr", seq)])
    FindOrfs(_args(fa, min_length=50, extract=True, protein=True)).run()
    out = capsys.readouterr().out
    lines = out.strip().split("\n")
    assert lines[0].startswith(">")
    # M + 99 A's, no stop ('*')
    assert lines[1] == "M" + "A" * 99
    assert "*" not in lines[1]


def test_find_orfs_json(tmp_path, capsys):
    seq = "AAAAAA" + CLEAN_ORF_FWD + "CCCCCC"
    fa = _write_fasta(tmp_path, [("chr", seq)])
    FindOrfs(_args(fa, min_length=50, format="json")).run()
    data = json.loads(capsys.readouterr().out)
    assert len(data) == 1
    assert data[0]["length_aa"] == 100
    assert data[0]["frame"] == "+1"


def test_find_orfs_no_orfs(tmp_path, capsys):
    # A short sequence with no ATGs
    fa = _write_fasta(tmp_path, [("chr", "CCCCCCCCCC")])
    FindOrfs(_args(fa, min_length=100)).run()
    out = capsys.readouterr().out.strip()
    assert out == "id\torf_id\tframe\tstart\tstop\tlength_nt\tlength_aa"


def test_find_orfs_rejects_zero_min_length(tmp_path):
    fa = _write_fasta(tmp_path, [("chr", "ATGAAA")])
    with pytest.raises(ValueError):
        FindOrfs(_args(fa, min_length=0))


def test_find_orfs_orf_without_stop(tmp_path, capsys):
    # ATG + 150 Ala codons with no stop codon; runs to end of sequence
    no_stop = "ATG" + "GCT" * 150
    fa = _write_fasta(tmp_path, [("chr", no_stop)])
    FindOrfs(_args(fa, min_length=100)).run()
    lines = capsys.readouterr().out.strip().split("\n")
    assert len(lines) == 2
    fields = lines[1].split("\t")
    # length_aa should be 151 (M + 150 Ala's), no stop codon consumed
    assert fields[6] == "151"
    # length_nt = 151 * 3 = 453
    assert fields[5] == "453"
