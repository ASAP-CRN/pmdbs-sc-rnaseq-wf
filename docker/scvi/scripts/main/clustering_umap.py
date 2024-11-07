import argparse
import scanpy
import leidenalg


parser = argparse.ArgumentParser(description="Annotate clusters")
parser.add_argument(
    "--adata-input", dest="adata_input", type=str, help="AnnData object for a dataset"
)
parser.add_argument(
    "--adata-output",
    dest="adata_output",
    type=str,
    help="Output file to save AnnData object to",
)
parser.add_argument(
    "--latent-key",
    dest="latent_key",
    type=str,
    default="X_scvi",
    help="Latent key to the scvi latent",
)

args = parser.parse_args()

### fixed parameters (TODO: make an argument)
n_neighbors = 15 #default
leiden_reslns = [0.05, 0.1, 0.2, 0.4]

# Set CPUs to use for parallel computing
scanpy._settings.ScanpyConfig.n_jobs = -1

adata = scanpy.read_h5ad(args.adata_input)  # type: ignore

# calculate neighbor graph on scVI latent
scanpy.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=args.latent_key)

# do leiden
for resolution in leiden_reslns:
    scanpy.tl.leiden(adata, resolution=resolution, key_added=f"leiden_res_{resolution:4.2f}")

scanpy.tl.umap(adata)

adata.write_h5ad(filename=args.adata_output, compression="gzip")
