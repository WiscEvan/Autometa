#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import os

import pandas as pd

from Bio import SeqIO

from autometa.common import coverage

assembly = "tests/data/metagenome.fna"


def test_from_spades_names():
    records = [rec for rec in SeqIO.parse(assembly, "fasta")]
    coverages = coverage.from_spades_names(records)
    assert len(coverages) == 3587
    assert coverages.index.name == "contig"


def test_make_length_table(tmpdir):
    temp_dir = tmpdir.mkdir("tempdir")
    len_table = coverage.make_length_table(
        fasta=assembly, out=os.path.join(temp_dir, "out")
    )
    assert type(len_table) is str
