"""
COPYRIGHT
Copyright 2020 Ian J. Miller, Evan R. Rees, Kyle Wolf, Siddharth Uppal,
Shaurya Chanana, Izaak Miller, Jason C. Kwan
This file is part of Autometa.
Autometa is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
Autometa is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.
You should have received a copy of the GNU Affero General Public License
along with Autometa. If not, see <http://www.gnu.org/licenses/>.
COPYRIGHT
Unit test for coverage.py
"""


import pytest
import pandas

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


@pytest.mark.skip
@pytest.mark.wip
def test_coverage_get_from_sam(small_metagenome, tmp_path):
    out = tmp_path / "covs_from_sam.tsv"
    df = coverage.get(fasta=small_metagenome, from_spades=False, out=out)
    assert df.index.name == "contig"
    assert "coverage" in df.columns
    assert out.exists()


@pytest.mark.skip
@pytest.mark.wip
def test_coverage_get_from_bam(small_metagenome, tmp_path):
    out = tmp_path / "covs_from_bam.tsv"
    df = coverage.get(fasta=small_metagenome, from_spades=False, out=out)
    assert df.index.name == "contig"
    assert "coverage" in df.columns
    assert out.exists()
