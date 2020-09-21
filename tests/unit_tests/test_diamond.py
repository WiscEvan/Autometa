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

Calculates coverage of contigs
"""

import pytest

from autometa.common.external import diamond

# TODO need to undertsnad more about attributes and use them here
def test_get_top_hit():
    hit1 = diamond.DiamondResult(
        qseqid="NODE_1379_length_8505_cov_228.207_1",
        sseqid="WP_046709068.1",
        pident=100,
        length=336,
        mismatch=0,
        gapopen=0,
        qstart=1,
        qend=336,
        sstart=1,
        send=336,
        evalue=2.7e-189,
        bitscore=671.4,
    )
    hit2 = diamond.DiamondResult(
        qseqid="NODE_1379_length_8505_cov_228.207_1",
        sseqid="WP_046918599.1",
        pident=100,
        length=336,
        mismatch=0,
        gapopen=0,
        qstart=1,
        qend=336,
        sstart=1,
        send=336,
        evalue=2.7e-189,
        bitscore=670.4,
    )
    hit3 = diamond.DiamondResult(
        qseqid="NODE_1379_length_8505_cov_228.207_1",
        sseqid="WP_037703909.1",
        pident=100,
        length=336,
        mismatch=0,
        gapopen=0,
        qstart=1,
        qend=336,
        sstart=1,
        send=336,
        evalue=2.7e-189,
        bitscore=675.4,
    )
    hit_list = [hit1, hit2, hit3]
    for hit in hit_list:
        if hit.sseqids == hit1.sseqids:
            hits = hit1
            continue
        hits = diamond.DiamondResult.__add__(hits, hit)
    # How to access bitscore?
    # hit1.sseqids.[sseqid].bitscore ?
    assert hits.get_top_hit() == "WP_037703909.1"


def test_parse():
    with pytest.raises(FileNotFoundError):
        diamond.parse(results="blast_results", bitscore_filter=0.9)
    # TODO mock diamond blast output file/mock blast call
    blast_results = (
        "/home/the_bio_informatician/Autometa/tests/data/dimaond_result_small"
    )
    with pytest.raises(ValueError):
        diamond.parse(results=blast_results, bitscore_filter=1.9)
    with pytest.raises(ValueError):
        diamond.parse(results=blast_results, bitscore_filter="TRUE")
    hits = diamond.parse(results=blast_results, bitscore_filter=0.9)
    assert type(hits) == dict


def test_add_taxids():
    blast_results = (
        "/home/the_bio_informatician/Autometa/tests/data/dimaond_result_small"
    )
    hits = diamond.parse(results=blast_results, bitscore_filter=0.9)
    with pytest.raises(FileNotFoundError):
        diamond.add_taxids(hits=hits, database="database")
