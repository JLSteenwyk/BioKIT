import json
from argparse import Namespace

import pytest

from biokit.services.text.kmer_frequency import KmerFrequency


def _args(fasta, kmer_size=3, canonical=False, verbose=False, format=None):
    return Namespace(
        fasta=fasta,
        kmer_size=kmer_size,
        canonical=canonical,
        verbose=verbose,
        format=format,
    )


def test_kmer_frequency_basic_tsv(tmp_path, capsys):
    fa = tmp_path / "seqs.fa"
    fa.write_text(">s\nAAAA\n")
    KmerFrequency(_args(str(fa), kmer_size=3)).run()
    lines = capsys.readouterr().out.strip().split("\n")
    # AAAA with k=3 → "AAA" at positions 0 and 1
    assert lines[0] == "kmer\tcount"
    assert lines[1] == "AAA\t2"


def test_kmer_frequency_multiple_kmers_sorted(tmp_path, capsys):
    fa = tmp_path / "seqs.fa"
    fa.write_text(">s\nATATATGCGC\n")
    KmerFrequency(_args(str(fa), kmer_size=2)).run()
    lines = capsys.readouterr().out.strip().split("\n")
    assert lines[0] == "kmer\tcount"
    rows = [line.split("\t") for line in lines[1:]]
    # Counts should be sorted descending
    counts = [int(r[1]) for r in rows]
    assert counts == sorted(counts, reverse=True)
    # AT appears 3x in ATATATGCGC (positions 0,2,4); TA 2x; GC 2x; TG 1x; CG 1x
    kmer_to_count = dict((r[0], int(r[1])) for r in rows)
    assert kmer_to_count["AT"] == 3
    assert kmer_to_count["TA"] == 2
    assert kmer_to_count["GC"] == 2
    assert kmer_to_count["TG"] == 1
    assert kmer_to_count["CG"] == 1


def test_kmer_frequency_case_insensitive(tmp_path, capsys):
    fa = tmp_path / "seqs.fa"
    fa.write_text(">s\nacgtacgt\n")
    KmerFrequency(_args(str(fa), kmer_size=4)).run()
    lines = capsys.readouterr().out.strip().split("\n")
    rows = [line.split("\t") for line in lines[1:]]
    kmer_to_count = dict((r[0], int(r[1])) for r in rows)
    assert kmer_to_count["ACGT"] == 2


def test_kmer_frequency_skips_ambiguous(tmp_path, capsys):
    fa = tmp_path / "seqs.fa"
    fa.write_text(">s\nACGNACG\n")
    KmerFrequency(_args(str(fa), kmer_size=3)).run()
    lines = capsys.readouterr().out.strip().split("\n")
    rows = [line.split("\t") for line in lines[1:]]
    kmers = [r[0] for r in rows]
    # Only ACG appears twice; CGN, GNA, NAC all contain N and are skipped
    assert kmers == ["ACG"]
    assert int(rows[0][1]) == 2


def test_kmer_frequency_canonical(tmp_path, capsys):
    fa = tmp_path / "seqs.fa"
    # ATGCAT has kmers: ATG, TGC, GCA, CAT (k=3)
    # ATG rc = CAT → canonical = min("ATG", "CAT") = "ATG"
    # TGC rc = GCA → canonical = min("TGC", "GCA") = "GCA"
    fa.write_text(">s\nATGCAT\n")
    KmerFrequency(_args(str(fa), kmer_size=3, canonical=True)).run()
    lines = capsys.readouterr().out.strip().split("\n")
    rows = [line.split("\t") for line in lines[1:]]
    kmer_to_count = dict((r[0], int(r[1])) for r in rows)
    # ATG + CAT collapse to ATG (count 2); TGC + GCA collapse to GCA (count 2)
    assert kmer_to_count == {"ATG": 2, "GCA": 2}


def test_kmer_frequency_aggregates_across_sequences(tmp_path, capsys):
    fa = tmp_path / "seqs.fa"
    fa.write_text(">s1\nAAAA\n>s2\nAAAA\n")
    KmerFrequency(_args(str(fa), kmer_size=3)).run()
    lines = capsys.readouterr().out.strip().split("\n")
    # Each sequence contributes 2 "AAA"s → total 4
    assert lines[1] == "AAA\t4"


def test_kmer_frequency_verbose_per_sequence(tmp_path, capsys):
    fa = tmp_path / "seqs.fa"
    fa.write_text(">s1\nAAAA\n>s2\nCCCC\n")
    KmerFrequency(_args(str(fa), kmer_size=3, verbose=True)).run()
    lines = capsys.readouterr().out.strip().split("\n")
    assert lines[0] == "id\tkmer\tcount"
    assert "s1\tAAA\t2" in lines
    assert "s2\tCCC\t2" in lines


def test_kmer_frequency_json(tmp_path, capsys):
    fa = tmp_path / "seqs.fa"
    fa.write_text(">s\nAAAA\n")
    KmerFrequency(_args(str(fa), kmer_size=3, format="json")).run()
    data = json.loads(capsys.readouterr().out)
    assert data == [{"kmer": "AAA", "count": 2}]


def test_kmer_frequency_rejects_zero_k(tmp_path):
    fa = tmp_path / "seqs.fa"
    fa.write_text(">s\nACGT\n")
    with pytest.raises(ValueError):
        KmerFrequency(_args(str(fa), kmer_size=0))


def test_kmer_frequency_sequence_shorter_than_k(tmp_path, capsys):
    fa = tmp_path / "seqs.fa"
    fa.write_text(">s\nAC\n")
    KmerFrequency(_args(str(fa), kmer_size=5)).run()
    out = capsys.readouterr().out.strip()
    # Only the header row, no k-mer rows
    assert out == "kmer\tcount"
