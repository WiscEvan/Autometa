import pytest
import subprocess

import pandas as pd


from autometa.common import coverage

# test_coverage_data.json
# Should be test_data.json keyed by file being tested
# e.g. coverage_test_data = variables["coverage"]


@pytest.fixture(name="small_metagenome")
def fixture_metagenome(variables, tmp_path):
    kmer_test_data = variables["kmers"]
    records = kmer_test_data["small_metagenome"]
    outlines = ""
    for record, seq in records.items():
        outlines += f"{record}\n{seq}\n"
    fpath = tmp_path / "small_metagenome.fna"
    with open(fpath, "w") as fh:
        fh.write(outlines)
    return fpath.as_posix()


@pytest.fixture(name="sam_alignment")
def fixture_alignment_sam(variables, tmp_path):
    coverage_test_data = variables["coverage"]
    alignment_records = coverage_test_data["sam"]
    outlines = ""
    for unique_id, values in alignment_records.items():
        values = "\t".join(values)
        outlines += f"{unique_id}\t{values}\n"
    fpath = tmp_path / "records.sam"
    with open(fpath, "w") as fh:
        fh.write(outlines)
    return fpath.as_posix()


@pytest.fixture(name="bam_alignment")
def fixture_alignment_bam(variables, sam_alignment, tmp_path):
    bam_fpath = tmp_path / "records.bam"
    cmd = f"samtools view -bS {sam_alignment} > {bam_fpath}"
    subprocess.run(
        cmd,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        shell=True,
        check=True,
    )
    return bam_fpath.as_posix()


@pytest.fixture(name="bed_alignment")
def fixture_alignment_bed(variables, tmp_path):
    coverage_test_data = variables["coverage"]
    alignment_records = coverage_test_data["bed"]
    outlines = ""
    for unique_id, values in alignment_records.items():
        unique_id = unique_id.rsplit("_", 1)[0]
        values = "\t".join(values)
        outlines += f"{unique_id}\t{values}\n"
    fpath = tmp_path / "records.bed"
    with open(fpath, "w") as fh:
        fh.write(outlines)
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


# TODO: Create fixtures for each set of data.
# Then group fixtures into a group and indirectly parametrize group
# i.e. alignments.sam, alignments.bam, alignments.bed, fwd_reads.fq, rev_reads.fq
# The external tools could be run or we could monkeypatch these.


def test_coverage_get_from_spades(small_metagenome, tmp_path):
    out = tmp_path / "covs_from_spades.tsv"
    df = coverage.get(fasta=small_metagenome, from_spades=True, out=out)
    assert df.index.name == "contig"
    assert "coverage" in df.columns
    assert out.exists()


def test_coverage_get_from_bed(small_metagenome, bed_alignment, tmp_path):
    out = tmp_path / "covs_from_bed.tsv"
    df = coverage.get(
        fasta=small_metagenome, from_spades=False, out=out, bed=bed_alignment
    )
    assert df.index.name == "contig"
    assert "coverage" in df.columns
    assert out.exists()


def test_coverage_get_from_sam(small_metagenome, sam_alignment, tmp_path):
    out = tmp_path / "covs_from_sam.tsv"
    df = coverage.get(
        fasta=small_metagenome, from_spades=False, out=out, sam=sam_alignment
    )
    assert df.index.name == "contig"
    assert "coverage" in df.columns
    assert out.exists()


def test_coverage_get_from_bam(small_metagenome, bam_alignment, tmp_path):
    out = tmp_path / "covs_from_bam.tsv"
    df = coverage.get(
        fasta=small_metagenome, from_spades=False, out=out, bam=bam_alignment
    )
    assert df.index.name == "contig"
    assert "coverage" in df.columns
    assert out.exists()


def test_get_ValueError(small_metagenome, tmp_path):
    out = tmp_path / "covs.tsv"
    with pytest.raises(ValueError):
        coverage.get(fasta=small_metagenome, from_spades=False, out=out)


def test_embed_df_already_exists(small_metagenome, df_exists_fpath, bed_alignment):
    coverage.get(fasta=small_metagenome, out=df_exists_fpath, bed=bed_alignment)
