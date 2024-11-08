import argparse
import scanpy as sc
import sys
import pandas as pd
from anndata import AnnData

sys.path.append("/opt/scripts/utility")
from helpers import score_cell_cycle, update_validation_metrics

def process_adata(adata: AnnData, markers:pd.DataFrame, batch_key:str) -> AnnData:
    """
    do feature selection and add PCA
    """
    
    ### fixed parameters (TODO: make an argument)
    n_components = 30
    target_sum = 1e4
    n_top_genes_all = 20_000


    # does this work with sparse uint8?
    adata.layers["counts"] = adata.X.copy()  # type: ignore

    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)

    # get cell cycle scores
    score_cell_cycle(adata, organism="human")

    ### MAKE SURE MARKER GENES ARE KEPT
    # defensive
    markers = markers[~markers.index.duplicated(keep="first")].rename_axis(index=None)

    batch_key = args.batch_key
    # WARNING: using 'sample' can cause loess to fail in the highly_variable_genes function
    # HACK: using 'batch_id' instead of 'sample' for now
    if batch_key == "sample":
        print("WARNING: using 'batch_id' instead of 'sample' for now")
        batch_key = "batch_id"

    # prefer pearson_residuals to seurat_v3
    # hvgs_full = scanpy.pp.highly_variable_genes(
    #     adata,
    #     n_top_genes=n_top_genes_all,
    #     batch_key=batch_key,
    #     flavor="seurat_v3",
    #     check_values=True,
    #     layer="counts",
    #     subset=False,
    #     inplace=False,
    # )
    hvgs_full = sc.experimental.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes_all,
        batch_key=batch_key,
        flavor="pearson_residuals",
        check_values=True,
        layer="counts",
        subset=False,
        inplace=False,
    )

    # hack to make sure we keep the marker genes
    hvgs_full.loc[markers.index, "highly_variable_nbatches"] = hvgs_full["highly_variable_nbatches"].max() + 1.0
    # Sort genes by how often they selected as hvg within each batch and
    # break ties with median rank of residual variance across batches
    hvgs_full.sort_values(
        ["highly_variable_nbatches", "highly_variable_rank"],
        ascending=[False, True],
        na_position="last",
        inplace=True,
    )
    hvgs_full = hvgs_full.iloc[: args.n_top_genes].index.to_list()
    adata = adata[:, adata.var.index.isin(hvgs_full)]

    # add PCA of log1p normalized
    sc.pp.pca(adata, n_comps=n_components)
    ## consider this alternate pca and avoid the normalization log1p?
    # scanpy.experimental.pp.normalize_pearson_residuals_pca(adata, *, theta=100, clip=None, n_comps=50, random_state=0, kwargs_pca=mappingproxy({}), mask_var=_empty, use_highly_variable=None, check_values=True, inplace=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Normalize seurat objects")

    parser.add_argument(
        "--adata-input", dest="adata_input", type=str, help="AnnData object for a dataset"
    )
    parser.add_argument(
        "--batch-key",
        dest="batch_key",
        type=str,
        help="Key in AnnData object for batch information",
    )
    parser.add_argument(
        "--adata-output",
        dest="adata_output",
        type=str,
        help="Output file to save AnnData object to",
    )
    parser.add_argument(
        "--n-top-genes",
        dest="n_top_genes",
        type=int,
        help="number of HVG genes to keep",
        default=3000,
    )
    parser.add_argument(
        "--marker-genes",
        dest="marker_genes",
        type=str,
        default="resources/celltype_marker_table.csv", # this doesn't seem accurate
        help="Path to marker_genes .csv file",
    )
    parser.add_argument(
        "--output-validation-file",
        dest="output_validation_file",
        type=str,
        help="Output file to write validation metrics to",
    )


    args = parser.parse_args()

    # Set CPUs to use for parallel computing
    sc._settings.ScanpyConfig.n_jobs = -1

    # 0. load data
    adata = sc.read_h5ad(args.adata_input)  # type: ignore
    # 1. load marker_genes
    # alternative way to get markers:
    # https://github.com/NIH-CARD/brain-taxonomy/blob/main/markers/cellassign_card_markers.csv
    markers = pd.read_csv(args.marker_genes, index_col=0)

    # 2. process data
    adata = process_adata(adata, markers, args.batch_key)
    
    # 3. save the filtered adata
    # save the filtered adata
    adata.write_h5ad(filename=args.adata_output, compression="gzip")

    # 4. update the validation metrics
    #######  validation metrics
    val_metrics = pd.read_csv(args.output_validation_file, index_col=0)
    output_metrics = update_validation_metrics(adata, "filter", val_metrics)
    # log the validation metrics
    output_metrics.to_csv(args.output_validation_file, index=True)

