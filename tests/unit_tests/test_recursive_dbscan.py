import pytest
import pandas as pd


from autometa.binning import recursive_dbscan
from autometa.common import kmers
from autometa.common.exceptions import TableFormatError


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
def fixture_embed_df(norm_df, tmpdir):
    embed_table = tmpdir.join("embed_table")
    df = kmers.embed(kmers=norm_df, out=embed_table)
    return df


@pytest.fixture(name="df_without_markers")
def fixture_empty_df():
    invalid_dict = {
        "contig": ["invalid_contig_1", "invalid_contig_2", "invalid_contig_3"],
        "markers": ["invalid_marker1", "invalid_marker2", "invalid_marker3"],
    }
    df = pd.DataFrame(invalid_dict)
    return df


# @pytest.mark.parametrize("method, verbose", [("dbscan", True), ("hdbscan", False)])
@pytest.mark.parametrize("verbose", [True, False])
# @pytest.mark.parametrize("verbose, reverse_ranks", [(True, True), (False, False)])
@pytest.mark.parametrize("method", ["dbscan", "hdbscan"])
def test_binning_bacteria(tmp_path, bacteria_markers, embed_kmers, method, verbose):
    out = tmp_path / "output_table.tsv"
    out = recursive_dbscan.binning(
        master=embed_kmers,
        markers=bacteria_markers,
        domain="bacteria",
        method=method,
        taxonomy=False,
        reverse_ranks=True,
        verbose=verbose,
    )
    assert len(out) == 50
    assert out.index.name == "contig"
    assert isinstance(out, pd.DataFrame)


@pytest.mark.wip
# @pytest.mark.parametrize("method, verbose", [("dbscan", True), ("hdbscan", False)])
# @pytest.mark.parametrize("verbose", [True, False])
# @pytest.mark.parametrize("verbose, reverse_ranks", [(True, True), (False, False)])
@pytest.mark.parametrize("method", ["dbscan", "hdbscan"])
def test_binning_taxonomy_true(tmp_path, bacteria_markers, embed_kmers, method):
    out = tmp_path / "output_table.tsv"
    out = recursive_dbscan.binning(
        master=embed_kmers,
        markers=bacteria_markers,
        domain="bacteria",
        method=method,
        taxonomy=True,
        reverse_ranks=True,
    )
    assert len(out) == 50
    assert out.index.name == "contig"
    assert isinstance(out, pd.DataFrame)


def test_binning_empty_markers_table(
    df_without_markers, embed_kmers,
):
    with pytest.raises(TableFormatError):
        recursive_dbscan.binning(
            master=embed_kmers,
            markers=df_without_markers,
            domain="bacteria",
            method="hdbscan",
            taxonomy=False,
        )


# @pytest.mark.parametrize("method", ["dbscan", "hdbscan"])
# # @pytest.mark.parametrize("verbose", [True, False])
# def test_binning_archaea(tmp_path, archaea_markers, embed_kmers, method):
#     out = tmp_path / "output_table.tsv"
#     out = recursive_dbscan.binning(
#         master=embed_kmers,
#         markers=archaea_markers,
#         domain="archaea",
#         method=method,
#         taxonomy=False,
#         verbose=False,
#     )
#     assert len(out) == 50
#     assert out.index.name == "contig"
#     # print(type(out))
#     assert isinstance(out, pd.DataFrame)


# @pytest.mark.skip
# @pytest.mark.parametrize(
#     "domain, verbose,reverse_ranks",
#     [("bacteria", True, True, pytest.param("bacteria_markers")), ("archaea", False, False, pytest.param("archaea_markers"))],
# )
# # @pytest.mark.parametrize(
# #     "markers_df",
# #     [pytest.param("bacteria_markers"), pytest.param("archaea_markers"),],
# #     indirect=True,
# # )
# @pytest.mark.parametrize(
#     "markers_df",
#     ["bacteria_markers", "archaea_markers"],
#     indirect=True,
# )
# # @pytest.mark.parametrize(
# #     "markers_df,domains",
# #     [pytest.param("bacteria_markers"), pytest.param("archaea_markers"),],
# #     indirect=True,
# # )
# def test_binning(markers_df, domain, verbose, reverse_ranks, embed_kmers, tmp_path):
#     out = tmp_path / "output_table.tsv"
#     out = recursive_dbscan.binning(
#         master=embed_kmers,
#         markers=markers_df,
#         domain=domain,
#         taxonomy=False,
#         verbose=verbose,
#         reverse_ranks=reverse_ranks,
#     )
#     assert len(out) == 50
#     assert out.index.name == "contig"
#     # print(type(out))
#     assert isinstance(out, pd.DataFrame)
