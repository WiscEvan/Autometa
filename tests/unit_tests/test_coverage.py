#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
Script containing wrapper functions for bedtools.
"""


import os
import pytest

import pandas as pd

from Bio import SeqIO
from unittest.mock import patch, MagicMock

from autometa.common import coverage

assembly = os.path.join("tests", "data", "metagenome.fna")
coverage_out = os.path.join("tests", "data", "coverage.tsv")


@pytest.fixture(scope="session")
def make_temp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("data")


def test_from_spades_names():
    records = [rec for rec in SeqIO.parse(assembly, "fasta")]
    coverages = coverage.from_spades_names(records)
    assert len(coverages) == 3587
    assert coverages.index.name == "contig"
    assert coverages.name == "coverage"
    assert type(coverages) is pd.core.series.Series
    assert coverages.dtypes == "float64"
    for cov in coverages.index:
        assert float(cov.split("_cov_")[-1]) == coverages[cov]


def test_make_length_table(make_temp_dir):
    # str() is needed as by default tmpdir_factory return -> "class py._path.local.LocalPath "
    # See https://stackoverflow.com/a/27034002/12671809
    out_fpath = str(make_temp_dir.join("out"))
    len_table = coverage.make_length_table(fasta=assembly, out=out_fpath)
    assert type(len_table) is str


@patch("os.path.getsize", return_value=2 * 1024 * 1024)
@patch("autometa.common.coverage.bowtie.align")
@patch("autometa.common.coverage.bowtie.build")
@patch("autometa.common.coverage.samtools.sort")
@patch("autometa.common.coverage.bedtools.genomecov")
@patch("autometa.common.coverage.bedtools.parse")
def test_get(
    patched_bedtools_parse,
    patched_bedtools_genomecov,
    patched_samtools_sort,
    patched_bowtie_build,
    patched_bowtie_align,
    patched_file_size,
    make_temp_dir,
):
    temp_out_fpath = make_temp_dir.join("cov_out")
    # Mocks of the files are needed to make sure that a Mock files "exists", when the test goes through coverage.py. Just putting a string in coverage.get() won't create a file and the test will fail while parsing through coverage.py

    bed_fpath = MagicMock(
        return_value="path_to_bed_file", name="path to input bed file",
    )
    bam_fpath = MagicMock(
        return_value="path_to_bam_file", name="path to input bam file",
    )
    sam_fpath = MagicMock(
        return_value="path_to_sam_file", name="path to input sam file",
    )
    fwd_reads_fpath = MagicMock(
        return_value="path_to_fwd_reads_file",
        name="path to fastq file having forward reads",
    )
    rev_reads_fpath = MagicMock(
        return_value="path_to_rev_reads_file",
        name="path to fastq file having reverse reads",
    )
    se_reads_fpath = MagicMock(
        return_value="path_to_se_reads_file",
        name="path to fastq file having single-end reads",
    )

    out_df = coverage.get(fasta=assembly, out=coverage_out)
    coverage.get(fasta=assembly, out=temp_out_fpath, bed=bed_fpath)
    coverage.get(fasta=assembly, out=temp_out_fpath, bam=bam_fpath)
    coverage.get(fasta=assembly, out=temp_out_fpath, sam=sam_fpath)
    coverage.get(
        fasta=assembly,
        out=temp_out_fpath,
        fwd_reads=fwd_reads_fpath,
        rev_reads=rev_reads_fpath,
    )
    coverage.get(
        fasta=assembly, out=temp_out_fpath, se_reads=se_reads_fpath,
    )

    assert len(out_df) == 3425
    assert out_df.index.name == "contig"
    assert type(out_df) is pd.DataFrame
    assert out_df.coverage.dtypes == "float64"
    assert list(out_df)[0] == "coverage"
    assert out_df.shape == (3425, 1)

    assert patched_bedtools_parse.call_count == 5
    assert patched_bedtools_genomecov.call_count == 4
    assert patched_samtools_sort.call_count == 3
    assert patched_bowtie_build.call_count == 2
    assert patched_bowtie_align.call_count == 2

    assert patched_file_size.call_count == 4
