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

Count, normalize and embed k-mers given nucleotide sequences
"""

import os
import pandas as pd
import pytest

from unittest.mock import patch

from autometa.common.metagenome import Metagenome
from autometa.common.external import prodigal

assembly = os.path.join("tests", "data", "metagenome.fna")


@pytest.fixture(scope="session")
def init_Metagenome(tmpdir_factory):
    temp_dir = tmpdir_factory.mktemp("data")
    mg = Metagenome(
        assembly=assembly,
        outdir=temp_dir,
        nucl_orfs_fpath="orfs.fna",
        prot_orfs_fpath="orfs.faa",
        fwd_reads=None,
        rev_reads=None,
        se_reads=None,
    )
    return mg


def test_Metagenome(init_Metagenome):
    assert init_Metagenome.largest_seq == "NODE_1_length_1389215_cov_225.275"
    assert init_Metagenome.size == 79468311
    assert init_Metagenome.nseqs == 3587
    assert round(init_Metagenome.length_weighted_gc, 2) == 54.34


def test_length_filter(init_Metagenome):
    out_fpath = os.path.join(init_Metagenome.outdir, "metagenome.filtered.fna")
    mg = init_Metagenome.length_filter(out=out_fpath, cutoff=3000)
    assert mg.largest_seq == "NODE_1_length_1389215_cov_225.275"
    assert mg.nseqs == 2059
    assert mg.size == 78034179
    assert round(mg.length_weighted_gc, 2) == 54.39

    with pytest.raises(ValueError):
        init_Metagenome.length_filter(out=out_fpath, cutoff="invalid_cutoff")
    with pytest.raises(ValueError):
        init_Metagenome.length_filter(out=out_fpath, cutoff=-6)
    with pytest.raises(FileExistsError):
        init_Metagenome.length_filter(out=out_fpath, cutoff="3000", force=False)


@patch(
    "autometa.common.external.prodigal.run",
    return_value=("nucls_out_fpath", "prots_out_fpath"),
)
def test_call_orfs(patched_progial_run, init_Metagenome):
    with pytest.raises(TypeError):
        init_Metagenome.call_orfs(force="invalid_input")
    with pytest.raises(TypeError):
        init_Metagenome.call_orfs(parallel="invalid_input")
    with pytest.raises(TypeError):
        init_Metagenome.call_orfs(cpus="invalid_input")

    nucls_fp, prots_fp = init_Metagenome.call_orfs()
    assert patched_progial_run.call_count == 1
    assert type(nucls_fp) is str
    assert type(prots_fp) is str
