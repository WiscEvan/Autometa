#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
COPYRIGHT
Copyright 2021 Ian J. Miller, Evan R. Rees, Kyle Wolf, Siddharth Uppal,
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

Autometa clustering evaluation benchmarking.

Script to benchmark Autometa clustering results using clustering evaluation metrics.
# Setting seed? See: https://stackoverflow.com/a/5837352/13118765
"""


import logging

import os
from typing import Dict, Iterable, NamedTuple, Union
import pandas as pd
from collections import namedtuple

from sklearn import metrics

from autometa.taxonomy.ncbi import NCBI

logger = logging.getLogger(__name__)

Labels = namedtuple("Labels", ["true", "pred"])
Targets = namedtuple("Targets", ["true", "pred", "target_names"])


def get_categorical_labels(
    predictions: str, reference: Union[str, pd.DataFrame]
) -> NamedTuple:
    pred_df = pd.read_csv(
        predictions, sep="\t", index_col="contig", usecols=["contig", "cluster"]
    )
    if not isinstance(reference, pd.DataFrame) and isinstance(reference, str):
        ref_df = pd.read_csv(
            reference,
            sep="\t",
            index_col="contig",
            usecols=["contig", "reference_genome"],
        )
    elif not isinstance(reference, pd.DataFrame) and not isinstance(reference, str):
        raise ValueError(f"reference is an invalid argument type: {type(reference)}")
    else:
        ref_df = reference
    # Modification here retrieves only contigs in binning dataset. This is not giving us the "complete" score.
    ref_df = ref_df[ref_df.index.isin(pred_df.index)]
    # Remove any contigs that were assigned as misassembled
    ref_df = ref_df[ref_df.reference_genome != "misassembled"]
    # Set reference_genome as category type so we can easily retrieve categorical labels
    # NOTE: These do not need to be a one-to-one correspondence to the prediction dataframe's categorical labels
    # (The exact label value does not matter for these metrics as we are only looking at the groupings [not classification])
    ref_df.reference_genome = ref_df.reference_genome.astype("category")
    # Assign "unclustered" to NA and convert 'cluster' column to categorical type
    unclustered_idx = pred_df[pred_df.cluster == "unclustered"].index
    pred_df.loc[unclustered_idx, "cluster"] = pd.NA
    pred_df.cluster = pred_df.cluster.astype("category")
    # Merge reference_assignments and predictions
    main_df = pd.merge(pred_df, ref_df, how="left", left_index=True, right_index=True)
    if main_df.empty:
        raise ValueError(
            "The provided reference community and predictions do not match!"
        )
    # Retrieve categorical values for each set of labels (truth and predicted)
    labels_true = main_df.reference_genome.cat.codes.values
    labels_pred = main_df.cluster.cat.codes.values
    return Labels(true=labels_true, pred=labels_pred)


def compute_clustering_metrics(labels: NamedTuple) -> Dict[str, float]:
    """Calculate various clustering performance metrics listed below.

    Note
    ----

    Some of these clustering performance evaluation metrics adjust for chance. This is dicussed in more detail in
    the sklearn user guide. - `Adjustment for chance in clustering performance evaluation <https://scikit-learn.org/stable/auto_examples/cluster/plot_adjusted_for_chance_measures.html#sphx-glr-auto-examples-cluster-plot-adjusted-for-chance-measures-py>`_

    This analysis suggests the most robust and reliable metrics to use as the number of clusters
    increases are adjusted rand index and adjusted mutual info score.

    Metrics
    -------

    * `adjusted mutual info score <https://scikit-learn.org/stable/modules/clustering.html#mutual-information-based-scores>`_
    * `geometric normalized mutual info score <https://scikit-learn.org/stable/modules/generated/sklearn.metrics.normalized_mutual_info_score.html#sklearn-metrics-normalized-mutual-info-score>`_
    * `adjusted rand index <https://scikit-learn.org/stable/modules/clustering.html#rand-index>`_
    * `homogeneity <https://scikit-learn.org/stable/modules/generated/sklearn.metrics.homogeneity_score.html#sklearn.metrics.homogeneity_score>`_
    * `completeness <https://scikit-learn.org/stable/modules/generated/sklearn.metrics.completeness_score.html#sklearn.metrics.completeness_score>`_
    * `V-measure <https://scikit-learn.org/stable/modules/clustering.html#homogeneity-completeness-and-v-measure>`_
    * `fowlkes-mallows score <https://scikit-learn.org/stable/modules/clustering.html#fowlkes-mallows-scores>`_


    Parameters
    ----------
    predictions : str
        Path to Autometa binning results.

    reference : str|pd.DataFrame
        Path to or dataframe of known reference community assignments to compare against Autometa binning results.

    Returns
    -------
    Dict[str, float]
        computed clustering evaluation metrics keyed by their metric

    Raises
    -------
    ValueError
        The input arguments are not the correct type (pd.DataFrame or str)
    ValueError
        The provided reference community and predictions do not match!

    """
    # Calculate cluster variety of cluster evaluation metrics with labels
    ari = metrics.adjusted_rand_score(labels_true=labels.true, labels_pred=labels.pred)
    homogeneity_score = metrics.homogeneity_score(labels.true, labels.pred)
    completeness_score = metrics.completeness_score(labels.true, labels.pred)
    v_measure = metrics.v_measure_score(labels.true, labels.pred)
    ami = metrics.adjusted_mutual_info_score(
        labels_true=labels.true, labels_pred=labels.pred
    )
    gnmi = metrics.normalized_mutual_info_score(
        labels_true=labels.true,
        labels_pred=labels.pred,
        average_method="geometric",
    )
    fowlkes_mallows_score = metrics.fowlkes_mallows_score(
        labels_true=labels.true, labels_pred=labels.pred
    )
    return {
        "adjusted mutual info score": ami,
        "geometric normalized mutual info score": gnmi,
        "adjusted rand score": ari,
        "homogeneity score": homogeneity_score,
        "completeness score": completeness_score,
        "V-measure": v_measure,
        "fowlkes-mallows score": fowlkes_mallows_score,
    }


def evaluate_clustering(predictions: Iterable, reference: str) -> pd.DataFrame:
    reference = pd.read_csv(reference, sep="\t", usecols=["contig", "reference_genome"])
    all_metrics = []
    for prediction in predictions:
        labels = get_categorical_labels(predictions=prediction, reference=reference)
        metrics = compute_clustering_metrics(labels)
        metrics.update({"dataset": os.path.basename(predictions)})
        all_metrics.append(metrics)
    df = pd.DataFrame(all_metrics).set_index("dataset")
    return df


def get_target_labels(
    prediction: str, reference: Union[str, pd.DataFrame], ncbi: Union[str, NCBI]
) -> namedtuple:
    """[summary]

    Note
    ----

    The exact label value matters for these metrics as we are
    looking at the available target labels for classification (not clustering)

    Parameters
    ----------
    prediction : str
        [description]
    reference : Union[str, pd.DataFrame]
        [description]
    ncbi : Union[str, NCBI]
        [description]

    Returns
    -------
    namedtuple
        [description]

    Raises
    ------
    ValueError
        [description]
    ValueError
        [description]
    """
    pred_df = pd.read_csv(
        prediction, sep="\t", index_col="contig", usecols=["contig", "taxid"]
    ).convert_dtypes()
    if not isinstance(reference, pd.DataFrame) and isinstance(reference, str):
        ref_df = (
            pd.read_csv(
                reference,
                sep="\t",
                index_col="contig",
                usecols=["contig", "taxid"],
            )
            .dropna(axis="index")
            .convert_dtypes()
        )
    elif not isinstance(reference, pd.DataFrame) and not isinstance(reference, str):
        raise ValueError(f"reference is an invalid argument type: {type(reference)}")
    else:
        ref_df = reference
    # Merge reference_assignments and predictions
    main_df = pd.merge(
        pred_df,
        ref_df,
        how="inner",
        left_index=True,
        right_index=True,
        suffixes=("_pred", "_true"),
    )
    if main_df.empty:
        raise ValueError(
            "The provided reference community and predictions do not match!"
        )
    # Convert any old taxids to new taxids from merged.dmp
    ncbi = NCBI(ncbi) if isinstance(ncbi, str) else ncbi
    main_df.taxid_pred = main_df.taxid_pred.map(
        lambda tid: ncbi.convert_taxid_dtype(tid)
    )
    main_df.taxid_true = main_df.taxid_true.map(
        lambda tid: ncbi.convert_taxid_dtype(tid)
    )
    # Create binary encoded matrix for multi-label classification metrics
    # First join strings s.t. taxid|taxid|... to be used with pd.str.get_dummies(sep='|')
    main_df["true_lineage"] = main_df.taxid_true.map(
        lambda t: "|".join(
            str(l.get("taxid")) for l in ncbi.lineage(t, canonical=False)
        )
    )
    main_df["pred_lineage"] = main_df.taxid_pred.map(
        lambda t: "|".join(
            str(l.get("taxid")) for l in ncbi.lineage(t, canonical=False)
        )
    )
    # Now create our binary encoded matrices (NOTE: These are multi-label classification matrices)
    y_true = main_df.true_lineage.str.get_dummies()
    y_pred = main_df.pred_lineage.str.get_dummies()
    # Now we need to ensure our columns have one-to-one correspondence with both dataframes
    # Retrieve columns in y_true but not in y_pred
    y_true_cols = y_true.loc[:, ~y_true.columns.isin(y_pred.columns)].columns
    # Retrieve columns in y_pred but not in y_true
    y_pred_cols = y_pred.loc[:, ~y_pred.columns.isin(y_true.columns)].columns
    # Now add these columns with 0's to reflect their absence in the other respective dataframe
    for col in y_true_cols:
        y_pred[col] = 0
    for col in y_pred_cols:
        y_true[col] = 0
    # Now we need to ensure all column indices correspond to each other between dataframes
    all_cols = y_true.columns.tolist()
    y_pred = y_pred[all_cols]
    return Targets(true=y_true, pred=y_pred, target_names=all_cols)


def compute_classification_metrics(labels: namedtuple) -> dict:
    return metrics.classification_report(
        y_true=labels.true,
        y_pred=labels.pred,
        target_names=labels.target_names,
        output_dict=True,
        zero_division=0,
    )


def evaluate_classification(
    predictions: Iterable,
    reference: str,
    ncbi: Union[str, NCBI],
    keep_averages=["weighted avg", "samples avg"],
) -> pd.DataFrame:
    # Read in community reference assignments
    reference = (
        pd.read_csv(
            reference, sep="\t", usecols=["contig", "taxid"], index_col="contig"
        )
        # Convert the taxid dtype to int
        .convert_dtypes()
        # Drop any contigs missing taxid classification
        .dropna(axis="index")
    )
    # Instantiate NCBI so we can coordinate taxids
    ncbi = NCBI(ncbi) if isinstance(ncbi, str) else ncbi
    all_metrics = []
    all_reports = []
    # Compute metrics for all provided predictions
    for prediction in predictions:
        # convert and merge taxids of reference assignments and predictions
        labels = get_target_labels(
            prediction=prediction, reference=reference, ncbi=ncbi
        )
        # Compute metrics across all canonical ranks
        report = compute_classification_metrics(labels)
        report.update({"dataset": os.path.basename(prediction)})
        all_reports.append(report)
        averages = {k: v for k, v in report.items() if k in keep_averages}
        for average, scores in averages.items():
            metrics = scores
            metrics.update(
                {"average": average, "dataset": os.path.basename(prediction)}
            )
            all_metrics.append(metrics)
    df = pd.DataFrame(all_metrics).set_index("dataset")
    return df, all_reports


def main():
    import argparse
    import logging as logger

    logger.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logger.DEBUG,
    )
    parser = argparse.ArgumentParser(
        description="Benchmark classification or clustering against reference assignments for the provided simulated/synthetic community.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--benchmark",
        help="Type of benchmarking to perform",
        choices={"classification", "clustering"},
        required=True,
    )
    parser.add_argument(
        "--predictions",
        help="Path to Autometa predictions (May specify multiple if they all correspond to the same `--reference` community ",
        metavar="filepath",
        nargs="*",
    )
    parser.add_argument(
        "--reference",
        help="Path to community reference assignments",
        metavar="filepath",
        required=True,
    )
    parser.add_argument(
        "--output-wide",
        help="Path to write benchmarking evaluation metrics (each metric receives its own column) (Default: `benchmark_type`_benchmarks.tsv.gz",
        metavar="filepath",
        required=False,
    )
    parser.add_argument(
        "--output-long",
        help="Path to write clustering evaluation metrics (metrics are stacked into one 'metric' column)",
        metavar="filepath",
        required=False,
    )
    parser.add_argument(
        "--output-classification-reports",
        help="Path to write classification evaluation reports",
        metavar="dirpath",
        required=False,
    )
    parser.add_argument(
        "--ncbi",
        help="Path to NCBI databases directory (Required with --benchmark=classification)",
        metavar="dirpath",
        required=False,
    )
    args = parser.parse_args()
    if args.benchmark == "clustering":
        df = evaluate_clustering(predictions=args.predictions, reference=args.reference)
    else:
        if not args.ncbi:
            raise ValueError(f"--ncbi is required for the classification benchmark!")
        ncbi = NCBI(args.ncbi)
        df, reports = evaluate_classification(
            predictions=args.predictions,
            reference=args.reference,
            ncbi=ncbi,
        )
        if args.output_classification_reports:
            if not os.path.isdir(
                args.output_classification_reports
            ) or not os.path.exists(args.output_classification_reports):
                os.makedirs(args.output_classification_reports)
                logger.info(
                    f"Created new directory: {args.output_classification_reports}"
                )
            for report in reports:
                dataset = report.pop("dataset")
                dataset = dataset.replace(".tsv", "").replace(".gz", "")
                dataset = f"{dataset}_classification_report.tsv.gz"
                report_filepath = os.path.join(
                    args.output_classification_reports, dataset
                )
                # Remove overall averages:
                avgs = [k for k in report if " avg" in k]
                for avg in avgs:
                    report.pop(avg)
                report_df = pd.DataFrame(report).transpose()
                report_df["name"] = report_df.index.map(lambda taxid: ncbi.name(taxid))
                report_df["rank"] = report_df.index.map(lambda taxid: ncbi.rank(taxid))
                report_df.index.name = "taxid"
                report_df.to_csv(report_filepath, sep="\t", index=True, header=True)
            logger.info(
                f"Wrote {len(reports):,} report(s) to {args.output_classification_reports}"
            )

    output_wide = (
        f"{args.benchmark}_benchmarks.tsv.gz"
        if not args.output_wide
        else args.output_wide
    )
    df.to_csv(output_wide, sep="\t", index=True, header=True)
    logger.info(f"Wrote {df.shape[0]} datasets metrics to {output_wide}")
    if args.output_long and not args.benchmark == "classification":
        # Write out stacked dataframe for visualization with `plot-cluster-evaluation-metrics.R`
        dff = df.stack()
        dff.index.name = ("dataset", "metric")
        dff.name = "score"
        dff = dff.to_frame().reset_index(level=1).rename(columns={"level_1": "metric"})
        dff.to_csv(args.output_long, sep="\t", index=True, header=True)
        logger.info(
            f"Wrote {dff.index.nunique()} datasets (stacked) metrics to {args.output_long}"
        )


if __name__ == "__main__":
    main()