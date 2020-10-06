import pytest
import subprocess
import os
import argparse

import pandas as pd


from autometa.common import coverage
from autometa.common.external import samtools


@pytest.fixture(name="metagenome")
def fixture_metagenome(variables, tmp_path):
    metagenome_test_data = variables["metagenome"]
    records = metagenome_test_data["assembly"]
    outlines = ""
    for record, seq in records.items():
        outlines += f"{record}\n{seq}\n"
    fpath = tmp_path / "metagenome.fna"
    with open(fpath, "w") as fh:
        fh.write(outlines)
    return fpath.as_posix()


@pytest.fixture(name="sam_alignment")
def fixture_alignment_sam(variables, tmp_path):
    coverage_test_data = variables["coverage"]
    outlines = coverage_test_data["sam"]
    fpath = tmp_path / "records.sam"
    with open(fpath, "w") as fh:
        fh.write(outlines)
    return fpath.as_posix()


@pytest.fixture(name="bam_alignment")
def fixture_alignment_bam(variables, sam_alignment, tmp_path):
    bam_fpath = tmp_path / "records.bam"
    samtools.sort(sam=sam_alignment, bam=bam_fpath)
    return bam_fpath.as_posix()


@pytest.fixture(name="bed_alignment")
def fixture_alignment_bed(variables, tmp_path):
    coverage_test_data = variables["coverage"]
    alignment_records = pd.read_json(coverage_test_data["bed"])
    fpath = tmp_path / "records.bed"
    alignment_records.to_csv(fpath, sep="\t", header=True, index=False)
    return fpath.as_posix()


@pytest.fixture(name="df_exists_fpath")
def fixture_df_without_contig_index_(tmp_path):
    df_dict = {
        "contig": ["contig_1", "contig_2", "contig_3"],
        "coverage": [1, 2, 3],
    }
    df = pd.DataFrame(df_dict)
    df_fpath = tmp_path / "invalid_df.tsv"
    df.to_csv(df_fpath, sep="\t")
    return df_fpath.as_posix()


@pytest.fixture(name="fwd_reads")
def fixture_fwd_reads(variables, tmp_path):
    coverage_test_data = variables["coverage"]
    outlines = coverage_test_data["fwd_reads"]
    fpath = tmp_path / "fwd_reads.fastq"
    with open(fpath, "w") as fh:
        fh.write(outlines)
    return fpath.as_posix()


@pytest.fixture(name="rev_reads")
def fixture_rev_reads(variables, tmp_path):
    coverage_test_data = variables["coverage"]
    outlines = coverage_test_data["rev_reads"]
    fpath = tmp_path / "rev_reads.fastq"
    with open(fpath, "w") as fh:
        fh.write(outlines)
    return fpath.as_posix()


def test_coverage_get_from_spades(metagenome, tmp_path):
    out = tmp_path / "covs_from_spades.tsv"
    df = coverage.get(fasta=metagenome, from_spades=True, out=out)
    assert df.index.name == "contig"
    assert "coverage" in df.columns
    assert out.exists()


def test_coverage_get_from_sam(metagenome, sam_alignment, tmp_path):
    out = tmp_path / "covs_from_sam.tsv"

    df = coverage.get(fasta=metagenome, from_spades=False, out=out, sam=sam_alignment)
    assert df.index.name == "contig"
    assert "coverage" in df.columns
    assert out.exists()


def test_coverage_get_from_bam(metagenome, bam_alignment, tmp_path):
    out = tmp_path / "covs_from_bam.tsv"
    df = coverage.get(fasta=metagenome, from_spades=False, out=out, bam=bam_alignment)
    assert df.index.name == "contig"
    assert "coverage" in df.columns
    assert out.exists()


def test_coverage_get_from_bed(metagenome, bed_alignment, tmp_path):
    out = tmp_path / "covs_from_bed.tsv"
    df = coverage.get(
        fasta=metagenome, from_spades=False, out=out.as_posix(), bed=bed_alignment
    )
    assert df.index.name == "contig"
    assert "coverage" in df.columns
    assert out.exists()


def test_coverage_get_from_reads(metagenome, fwd_reads, rev_reads, tmp_path):
    out = tmp_path / "covs_from_bed.tsv"
    df = coverage.get(
        fasta=metagenome,
        from_spades=False,
        out=out.as_posix(),
        fwd_reads=fwd_reads,
        rev_reads=rev_reads,
    )
    assert df.index.name == "contig"
    assert "coverage" in df.columns
    assert out.exists()


def test_get_ValueError(metagenome, tmp_path):
    out = tmp_path / "covs.tsv"
    with pytest.raises(ValueError):
        coverage.get(fasta=metagenome, from_spades=False, out=out)


def test_embed_df_already_exists(metagenome, df_exists_fpath, bed_alignment):
    coverage.get(fasta=metagenome, out=df_exists_fpath, bed=bed_alignment)


use_spades = [True, False]


@pytest.fixture(name="mock_parser", params=use_spades)
def fixture_mock_parser(
    metagenome,
    sam_alignment,
    bam_alignment,
    bed_alignment,
    fwd_reads,
    rev_reads,
    monkeypatch,
    tmp_path,
    request,
):
    def return_mock_parser(*args, **kwargs):
        return MockParser()

    class MockParseArgs:
        def __init__(
            self,
            metagenome,
            sam_alignment,
            bam_alignment,
            bed_alignment,
            fwd_reads,
            rev_reads,
            out,
        ):
            self.assembly = metagenome
            self.fwd_reads = fwd_reads
            self.rev_reads = rev_reads
            self.sam = sam_alignment
            self.bam = bam_alignment
            self.lengths = "lengths.tsv"
            self.bed = bed_alignment
            self.cpus = 2
            self.out = out
            self.from_spades = request.param

    class MockParser:
        def add_argument(self, *args, **kwargs):
            pass

        def parse_args(self):
            out = tmp_path / "binning.tsv"
            return MockParseArgs(
                metagenome,
                sam_alignment,
                bam_alignment,
                bed_alignment,
                fwd_reads,
                rev_reads,
                out,
            )

    # Defining the MockParser class to represent parser
    monkeypatch.setattr(argparse, "ArgumentParser", return_mock_parser, raising=True)


def test_coverage_main(monkeypatch, mock_parser):
    coverage.main()
