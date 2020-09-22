import pytest
import pandas as pd


from autometa.binning import recursive_dbscan
from autometa.common import kmers


@pytest.fixture(name="bacteria_markers")
def fixture_bacteria_markers(variables, tmp_path):
    markers_test_data = variables["markers"]
    df = pd.read_json(markers_test_data["bacteria"])
    df.set_index("contig", inplace=True)
    # out = tmp_path / "bacteria.markers.tsv"
    # df.to_csv(out, sep="\t", index=True, header=True)
    return df


@pytest.fixture(name="archaea_markers")
def fixture_archaea_markers(variables, tmp_path):
    markers_test_data = variables["markers"]
    df = pd.read_json(markers_test_data["archaea"])
    df.set_index("contig", inplace=True)
    # out = tmp_path / "archaea.markers.tsv"
    # df.to_csv(out, sep="\t", index=True, header=True)
    return df


@pytest.fixture(name="norm_df")
def fixture_norm_df(variables):
    kmer_test_data = variables["kmers"]
    df = pd.read_json(kmer_test_data["norm_df"])
    df.set_index("contig", inplace=True)
    return df


@pytest.fixture(name="embed_kmers")
def test_embed(norm_df, tmpdir):
    embed_table = tmpdir.join("embed_table")
    df = kmers.embed(kmers=norm_df, out=embed_table)
    return df


@pytest.fixture(name="markers_df")
def fixture_markers_fpath(request, bacteria_markers, archaea_markers):
    kingdom = request.param
    markers_df = {
        "bacteria_markers": bacteria_markers,
        "archaea_markers": archaea_markers,
    }
    return markers_df[kingdom]


# @pytest.mark.parametrize(
#     "markers_df, domain",
#     [("bacteria_markers", "bacteria"), ("archaea_markers", "archaea"),],
# )

# @pytest.mark.parametrize(
#     "markers_df", ["bacteria_markers", "archaea_markers"],
# )
# @pytest.mark.skip
@pytest.mark.parametrize("domain, verbose", [("bacteria", True), ("archaea", False)])
@pytest.mark.parametrize(
    "markers_df",
    [pytest.param("bacteria_markers"), pytest.param("archaea_markers"),],
    indirect=True,
)
def test_resursive_dbscan(markers_df, domain, verbose, embed_kmers, tmp_path):
    out_passed = tmp_path / "passed_contigs.tsv"
    out_failed = tmp_path / "failed_contigs.tsv"
    out_passed, out_failed = recursive_dbscan.recursive_dbscan(
        table=embed_kmers,
        markers_df=markers_df,
        domain=domain,
        completeness_cutoff=20,
        purity_cutoff=90,
        verbose=verbose,
    )
    assert len(out_failed) + len(out_passed) == 50
    assert out_passed.empty
    assert out_failed.index.name == "contig"


# @pytest.mark.skip
@pytest.mark.parametrize(
    "domain, verbose", [("bacteria", True), ("archaea", False)],
)
@pytest.mark.parametrize(
    "markers_df",
    [pytest.param("bacteria_markers"), pytest.param("archaea_markers"),],
    indirect=True,
)
def test_recursive_hdbscan(markers_df, domain, verbose, embed_kmers, tmp_path):
    out_passed = tmp_path / "passed_contigs.tsv"
    out_failed = tmp_path / "failed_contigs.tsv"
    out_passed, out_failed = recursive_dbscan.recursive_hdbscan(
        table=embed_kmers,
        markers_df=markers_df,
        domain="bacteria",
        completeness_cutoff=20,
        purity_cutoff=90,
        verbose=verbose,
    )
    assert len(out_failed) + len(out_passed) == 50
    assert out_passed.empty
    assert out_failed.index.name == "contig"


# @pytest.mark.skip
@pytest.mark.parametrize(
    "domain, verbose,reverse_ranks",
    [("bacteria", True, True), ("archaea", False, False)],
)
@pytest.mark.parametrize(
    "markers_df",
    [pytest.param("bacteria_markers"), pytest.param("archaea_markers"),],
    indirect=True,
)
def test_binning(markers_df, domain, verbose, reverse_ranks, embed_kmers, tmp_path):
    out = tmp_path / "output_table.tsv"
    out = recursive_dbscan.binning(
        master=embed_kmers,
        markers=markers_df,
        domain=domain,
        taxonomy=False,
        verbose=verbose,
        reverse_ranks=reverse_ranks,
    )
