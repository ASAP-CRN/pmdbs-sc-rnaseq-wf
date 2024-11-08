import argparse
import scanpy as sc
from anndata import concat as ad_concat
from anndata import AnnData
import pandas as pd
import sys

sys.path.append("/opt/scripts/utility")
from helpers import get_validation_metrics


def merge_adata(adata_objects_fofn: str) -> AnnData:
    adatas = {}
    #  note that the sample id should be the official ASAP_samples
    with open(adata_objects_fofn, "r") as file:
        for sample in file:
            # Check that sample line is not empty
            if sample.strip():
                columns = sample.strip().split("\t")
                sample_name, file_path = columns
                raw = sc.read_h5ad(file_path)
                adatas[sample_name] = raw
    # we could subset to the top_genes here before concat if we have memory issues (e.g. whole dataset harmonization.)
    adata = ad_concat(merge="same", uns_merge="same", index_unique="_", adatas=adatas)

    return adata


def gen_qc_plots(adata: AnnData):
    # Fixed parameters
    metrics = [
        "n_genes_by_counts",
        "total_counts",
        "pct_counts_mt",
        "pct_counts_rb",
        "doublet_score",
    ]

    # Set CPUs to use for parallel computing
    sc._settings.ScanpyConfig.n_jobs = -1

    sc.settings.verbosity = 1
    sc.settings.figdir = "plots/"
    sc.settings.set_figure_params(
        dpi=100, fontsize=10, dpi_save=300, format="png", figsize=("12", "8")
    )  # type: ignore

    sc.settings.verbosity = 1
    sc.settings.figdir = "plots/"
    sc.settings.set_figure_params(
        dpi=100, fontsize=10, dpi_save=300, format="png", figsize=("12", "8")
    )  # type: ignore

    for metric in metrics:  # type: ignore
        sc.pl.violin(adata, keys=metric, size=0, save="".join("_" + metric))


def main(args: argparse.Namespace):
    """
    basic logic with args as input

    """
    adata = merge_adata(args.adata_objects_fofn)
    gen_qc_plots(adata)

    # export concatenated data.
    adata.write_h5ad(filename=args.adata_output, compression="gzip")
    # save metadata
    adata.obs.to_csv(args.output_metadata_file, index=True)

    #######  validation metrics
    val_metrics = get_validation_metrics(adata, "concatenation")
    # log the validation metrics
    val_metrics.to_csv(args.output_validation_file, index=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Call doublets")
    parser.add_argument(
        "--adata-objects-fofn",
        dest="adata_objects_fofn",
        type=str,
        help="Newline-delimited TSV of sample names and paths to the set of input adata objects (file-of-filenames)",
    )
    parser.add_argument(
        "--adata-output",
        dest="adata_output",
        type=str,
        help="Output file to save AnnData object to",
    )
    parser.add_argument(
        "--output-metadata-file",
        dest="output_metadata_file",
        type=str,
        help="Output file to write metadata to",
    )
    parser.add_argument(
        "--output-validation-file",
        dest="output_validation_file",
        type=str,
        help="Output file to write validation metrics to",
    )

    args = parser.parse_args()
    main(args)
