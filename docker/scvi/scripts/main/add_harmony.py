import argparse
import scanpy

parser = argparse.ArgumentParser(description="Add PCA and Harmony integration")
parser.add_argument(
    "--batch-key",
    dest="batch_key",
    type=str,
    help="Key in AnnData object for batch information",
)
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
    "--output-metadata-file",
    dest="output_metadata_file",
    type=str,
    help="Output file to write metadata to",
)
args = parser.parse_args()
 
# Fixed parameters
n_comps = 30

# Set CPUs to use for parallel computing
scanpy._settings.ScanpyConfig.n_jobs = -1

adata = scanpy.read_h5ad(args.adata_input)  # type: ignore

# make sure we have PCA (should already be there)
if "X_pca" not in adata.obsm:
    scanpy.pp.pca(adata, n_comps=n_comps)

# add harmony
if "X_pca_harmony" not in adata.obsm:
    scanpy.external.pp.harmony_integrate(adata, args.batch_key)
    # 9. write_h5ad

adata.write_h5ad(filename=args.adata_output, compression="gzip")

# save metadata
metatable = adata.obs 
metatable['UMAP_1'] = adata.obsm['X_umap'][:,0]
metatable['UMAP_2'] = adata.obsm['X_umap'][:,1]

metatable.to_csv(args.output_metadata_file, index=True) 