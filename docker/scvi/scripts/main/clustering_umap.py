import argparse
import scanpy as sc
from anndata import AnnData


SCANVI_LATENT_KEY = "X_scANVI"
SCVI_LATENT_KEY = "X_scVI"

SCANVI_PREDICTIONS_KEY = "C_scANVI"

def get_cluster_umap(adata: AnnData, latent_key: str) -> AnnData:
    ### fixed parameters (TODO: make an argument)
    n_neighbors = 15  # default
    leiden_reslns = [0.05, 0.1, 0.2, 0.4]
    # Set CPUs to use for parallel computing
    sc._settings.ScanpyConfig.n_jobs = -1

    # calculate neighbor graph on scVI latent
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=latent_key)
    # do leiden
    # change to igraph for future default compatibility
    for resolution in leiden_reslns:
        sc.tl.leiden(
            adata, resolution=resolution, key_added=f"leiden_res_{resolution:4.2f}", flavor="igraph", n_iterations=2, directed = False
        )
    sc.tl.umap(adata)
    return adata


def main(args: argparse.Namespace):
    """
    basic logic with args as input

    """
    adata = sc.read_h5ad(args.adata_input)  # type: ignore
    # Hard code latent_key
    latent_key = SCVI_LATENT_KEY
    adata = get_cluster_umap(adata, latent_key)
    adata.write_h5ad(filename=args.adata_output, compression="gzip")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate clusters")
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

    args = parser.parse_args()
    main(args)
