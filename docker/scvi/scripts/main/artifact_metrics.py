import argparse
import scanpy as sc
import numpy as np
from scib_metrics.benchmark import Benchmarker, BioConservation
from scib_metrics.nearest_neighbors import NeighborsResults
import faiss
from pathlib import Path

from anndata import AnnData
from pandas import DataFrame


# TODO: add these functions to utility/helpers.py
def faiss_hnsw_nn(X: np.ndarray, k: int):
    """Gpu HNSW nearest neighbor search using faiss.

    See https://github.com/nmslib/hnswlib/blob/master/ALGO_PARAMS.md
    for index param details.
    """
    X = np.ascontiguousarray(X, dtype=np.float32)
    res = faiss.StandardGpuResources()
    M = 32
    index = faiss.IndexHNSWFlat(X.shape[1], M, faiss.METRIC_L2)
    gpu_index = faiss.index_cpu_to_gpu(res, 0, index)
    gpu_index.add(X)
    distances, indices = gpu_index.search(X, k)
    del index
    del gpu_index
    # distances are squared
    return NeighborsResults(indices=indices, distances=np.sqrt(distances))


def faiss_brute_force_nn(X: np.ndarray, k: int):
    """Gpu brute force nearest neighbor search using faiss."""
    X = np.ascontiguousarray(X, dtype=np.float32)
    res = faiss.StandardGpuResources()
    index = faiss.IndexFlatL2(X.shape[1])
    gpu_index = faiss.index_cpu_to_gpu(res, 0, index)
    gpu_index.add(X)
    distances, indices = gpu_index.search(X, k)
    del index
    del gpu_index
    # distances are squared
    return NeighborsResults(indices=indices, distances=np.sqrt(distances))


def get_artifact_metrics(
    adata: AnnData, batch_key: str, label_key: str, scib_report_dir: Path | str
) -> DataFrame:
    # Fixed parameters
    n_comps = 30

    # Set CPUs to use for parallel computing
    sc._settings.ScanpyConfig.n_jobs = -1

    # these should be there...
    if "X_pca" not in adata.obsm:
        sc.pp.pca(adata, n_comps=n_comps)

    if "X_pca_harmony" not in adata.obsm:
        sc.external.pp.harmony_integrate(adata, "sample")

    adata.obsm["Unintegrated"] = adata.obsm["X_pca"]
    biocons = BioConservation(isolated_labels=False)
    bm = Benchmarker(
        adata,
        batch_key=batch_key,
        label_key=label_key,
        embedding_obsm_keys=["Unintegrated", "X_scvi", "X_pca_harmony"],
        pre_integrated_embedding_obsm_key="X_pca",
        bio_conservation_metrics=biocons,
        n_jobs=-1,
    )
    bm.prepare(neighbor_computer=faiss_brute_force_nn)
    bm.benchmark()

    bm.plot_results_table(min_max_scale=False, save_dir=scib_report_dir)
    df = bm.get_results(min_max_scale=False)

    return df


def main(args: argparse.Namespace):
    """
    basic logic with args as input

    """
    adata = sc.read_h5ad(args.adata_input)  # type: ignore

    # save the results
    report_dir = Path.cwd() / args.scib_report_dir
    if not report_dir.exists():
        report_dir.mkdir(parents=True, exist_ok=True)

    report_df = get_artifact_metrics(adata, args.batch_key, args.label_key, report_dir)
    report_df.to_csv((report_dir / "scib_report.csv"), index=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute `scib` metrics on final artefacts"
    )
    parser.add_argument(
        "--label-key",
        dest="label_key",
        type=str,
        default="cell_type",
        help="key to reference our 'cell_type' labels",
    )
    parser.add_argument(
        "--batch-key",
        dest="batch_key",
        type=str,
        help="Key in AnnData object for batch information",
    )
    parser.add_argument(
        "--adata-input",
        dest="adata_input",
        type=str,
        help="AnnData object for a dataset",
    )
    parser.add_argument(
        "--output-report-dir",
        dest="scib_report_dir",
        type=str,
        help="Output folder to save `scib` report",
    )

    args = parser.parse_args()
    main(args)
