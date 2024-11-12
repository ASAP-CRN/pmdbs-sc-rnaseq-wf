# import muon.pp.filter_obs as filter_obs ???
import muon as mu
import scanpy as sc
import argparse
import pandas as pd
import sys
from anndata import AnnData

sys.path.append("/opt/scripts/utility")
from helpers import update_validation_metrics


def filter_adata(adata: AnnData) -> AnnData:
    # TODO: make these cutoffs arguments...
    ## Fixed Parameters
    pct_counts_mt_max = 10
    doublet_score_max = 0.2
    total_counts_limits = [100, 100000]
    n_genes_by_counts_limits = [100, 10000]

    # muon api is better than the scanpyt api for this...
    mu.pp.filter_obs(adata, "pct_counts_mt", lambda x: x <= pct_counts_mt_max)
    mu.pp.filter_obs(adata, "doublet_score", lambda x: x < doublet_score_max)
    mu.pp.filter_obs(
        adata,
        "total_counts",
        lambda x: (x >= total_counts_limits[0]) & (x <= total_counts_limits[1]),
    )
    mu.pp.filter_obs(
        adata,
        "n_genes_by_counts",
        lambda x: (x >= n_genes_by_counts_limits[0])
        & (x <= n_genes_by_counts_limits[1]),
    )

    return adata


def main(args: argparse.Namespace):
    """
    basic logic with args as input

    """
    # Set CPUs to use for parallel computing
    sc._settings.ScanpyConfig.n_jobs = -1

    # 0. load data
    adata = sc.read_h5ad(args.adata_input, backed="r")  # type: ignore
    # 1. filter data
    adata = filter_adata(adata)
    # 2. save the filtered adata
    #    save the filtered adata
    adata.write_h5ad(filename=args.adata_output, compression="gzip")
    # 3. update the validation metrics
    #    ######  validation metrics
    val_metrics = pd.read_csv(args.output_validation_file, index_col=0)
    output_metrics = update_validation_metrics(adata, "filter", val_metrics)
    # log the validation metrics
    output_metrics.to_csv(args.output_validation_file, index=True)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Filter")
    parser.add_argument(
        "--adata-input",
        dest="adata_input",
        type=str,
        help="AnnData object for a dataset",
    )
    parser.add_argument(
        "--adata-output",
        dest="adata_output",
        type=str,
        help="Output file to save AnnData object to",
    )
    # TODO: add filter parameters as arguments
    parser.add_argument(
        "--output-validation-file",
        dest="output_validation_file",
        type=str,
        help="Output file to write validation metrics to",
    )

    args = parser.parse_args()
    main(args)
