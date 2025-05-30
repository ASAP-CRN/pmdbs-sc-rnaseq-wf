# import muon.pp.filter_obs as filter_obs ???
import muon as mu
import scanpy as sc
import argparse
import pandas as pd
import sys
from anndata import AnnData

sys.path.append("/opt/scripts/utility")
from helpers import update_validation_metrics

from pathlib import Path
import os


# os.environ["AIBS_BKP_USE_TORCH"] = "false"

# # os.environ["NUMEXPR_NUM_THREADS"] = "1"
# # os.environ["MKL_NUM_THREADS"] = "1"
# # os.environ["OMP_NUM_THREADS"] = "1"

# HOME = Path.home()

# ROOT = HOME / "Projects/ASAP/cell_type_mapper"


# TMP_DIR = ROOT / "tmp"
# if not TMP_DIR.exists():
#     TMP_DIR.mkdir()



# CHUNK_SIZE = 40000
# N_RUNNERS_UP = 5
# RNG_SEED = 11235813
# N_PROCESSORS = 8
# MAX_GB = 48.0


# #
# FILE_ROOT = "asap-cohort.merged_adata_object"
# DATE = "20250129"


# REFERENCE = "SEAAD"

# EXTENDED_RESULTS = f"{FILE_ROOT}.mmc.{REFERENCE}.{DATE}.json"
# LOG_FILE = f"{FILE_ROOT}.mmc.{REFERENCE}_log.{DATE}.txt"
# CSV_RESULTS = f"{FILE_ROOT}.mmc.{REFERENCE}_results.{DATE}.csv"


# RESULTS_DIR = ROOT / f"MMC.{REFERENCE}_RESULTS"
# if not RESULTS_DIR.exists():
#     RESULTS_DIR.mkdir()

# PRECOMPUTED_STATS = (
#     f"examples/data/abc_atlas_data/precomputed_stats.20231120.sea_ad.MTG.h5"
# )

# load resuults.
# results header is the first 4 lines
def read_csv_results(csv_results_path: str | Path) -> pd.DataFrame:
    """
    Read the results file and return a pandas DataFrame
    """
    ## HEADER
    # metadata = asap-cohort.merged_adata_object.mmc.seaad_results.20250129.json
    # taxonomy hierarchy = ["CCN20230505_CLAS", "CCN20230505_SUBC", "CCN20230505_SUPT"]
    # readable taxonomy hierarchy = ["class", "subclass", "supertype"]
    # algorithm: 'hierarchical'; codebase: http://github.com/AllenInstitute/cell_type_mapper; version: v1.4.2
    # cell_id,class_label,class_name,class_bootstrapping_probability,class_correlation_coefficient,class_aggregate_probability,subclass_label,subclass_name,subclass_bootstrapping_probability,subclass_correlation_coefficient,subclas

    results = pd.read_csv(csv_results_path, header=4)
    results.set_index("cell_id", inplace=True)
    results.index.name = None
    return results



def summarize_mmc_results(mmc_results : pd.DataFrame):

    # First get "classes"
    #  define row indices for each class
    gabaergic = mmc_results["class_name"] == "Neuronal: GABAergic"
    glutamatergic = mmc_results["class_name"] == "Neuronal: Glutamatergic"
    non_neuronal = mmc_results["class_name"] == "Non-neuronal and Non-neural"

    mmc_results.loc[gabaergic, "phenotype"] = "GABAergic"
    mmc_results.loc[glutamatergic, "phenotype"] = "Glutamatergic"
    mmc_results.loc[non_neuronal, "phenotype"] = mmc_results.loc[non_neuronal, "subclass_name"]

    mmc_results.loc[glutamatergic, "rho"] = mmc_results.loc[glutamatergic,  "class_correlation_coefficient"]
    mmc_results.loc[gabaergic, "rho"] = mmc_results.loc[gabaergic, "class_correlation_coefficient"]
    mmc_results.loc[non_neuronal, "rho"] = mmc_results.loc[non_neuronal, "subclass_correlation_coefficient"]

    mmc_results.loc[glutamatergic, "prob"] = mmc_results.loc[glutamatergic, "class_bootstrapping_probability"]
    mmc_results.loc[gabaergic, "prob"] = mmc_results.loc[gabaergic, "class_bootstrapping_probability"]
    mmc_results.loc[non_neuronal, "prob"] = mmc_results.loc[non_neuronal, "subclass_bootstrapping_probability"]


    mmc_results["cell_type"] = mmc_results["phenotype"]

    # change the phenotype to unknown if the  correlation or bootstrap probability < 0.5
    mmc_results.loc[mmc_results["rho"] < 0.5, "cell_type"] = "unknown"
    mmc_results.loc[mmc_results["prob"] < 0.5, "cell_type"] = "unknown"


    return mmc_results[["cell_type", "phenotype", "rho", "prob","class_name", "subclass_name", "supertype_name" ]]



def main(args: argparse.Namespace):
    """
    basic logic with args as input

    """

    # load MMC results
    results = read_csv_results(args.mmc_results) 

    results = summarize_mmc_results(results)

    # save the results to parquet file.. or feather?
    results.to_parquet(args.cell_type_output)

    # load the adata
    adata = sc.read_h5ad(args.adata_input)

    adata.obs = adata.obs.merge(results, left_index=True, right_index=True)

    # save the adata
    adata.write_h5ad(filename=args.adata_output, compression="gzip")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Transcriptional Phenotype (MMC)")
    parser.add_argument(
        "--adata-input",
        dest="adata_input",
        type=str,
        help="AnnData object for a dataset",
    )
    parser.add_argument(
        "--output-cell-types-file",
        dest="cell_type_output",
        type=str,
        help="Output file to write celltypes to",
    )
    parser.add_argument(
        "--adata-output",
        dest="adata_output",
        type=str,
        help="Output file to save AnnData object to",
    )
    parser.add_argument(
        "--mmc-results",
        dest="mmc_results",
        type=str,
        help="Path to MMC result CSV",
    )
    
    args = parser.parse_args()
    main(args)
