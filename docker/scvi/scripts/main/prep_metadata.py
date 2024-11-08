import argparse
import sys
import scanpy as sc
from anndata import AnnData
sys.path.append("/opt/scripts/utility")
from helpers import anndata_from_h5

def prep_metadata(adata:AnnData, sample_id:str, batch_id:str, dataset_id:str, team_id:str):
    # Fixed parameters
    # Set CPUs to use for parallel computing
    sc._settings.ScanpyConfig.n_jobs = -1

    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["rb"] = adata.var_names.str.startswith(("RPL", "RPS"))
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["rb", "mt"], percent_top=None, log1p=False, inplace=True
    )

    # add doublet score
    sc.external.pp.scrublet(adata)
    # keep the metadata clean ad drop predicted_doublet... rely on the probability of dublet i.e. doublet_score
    adata.obs.drop("predicted_doublet", axis=1, inplace=True)

    # add metadata
    adata.obs["sample"] = sample_id
    adata.obs["batch"] = batch_id
    adata.obs["team"] = team_id
    adata.obs["dataset"] = dataset_id
    # adata.obs['batch_id'] = args.project+args.batch
    adata.obs["batch_id"] = f"{team_id}_{dataset_id}_{batch_id}"  #
    adata.write_h5ad(filename=args.adata_output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Preprocess")
    parser.add_argument(
        "--adata-input", dest="adata_input", type=str, help="AnnData object for a dataset"
    )
    # toplevel metadata to add to the adata
    parser.add_argument("--sample-id", dest="sample_id", type=str, help="Sample/dataset ID")
    parser.add_argument(
        "--batch",
        dest="batch",
        type=str,
        help="Batch from which the sample/dataset originated",
    )
    # "project" is ambiguous... change to dataset... add tea,
    parser.add_argument("--dataset", dest="dataset", type=str, help="Dataset ID")
    parser.add_argument("--team", dest="team", type=str, help="Team ID")
    parser.add_argument(
        "--adata-output",
        dest="adata_output",
        type=str,
        help="Output file to save AnnData object to",
    )

    args = parser.parse_args()

    # load the data from cellbender output
    adata = anndata_from_h5(args.adata_input)
    adata = prep_metadata(args.sample_id, args.batch, args.dataset, args.team)
    adata.write_h5ad(filename=args.adata_output)
