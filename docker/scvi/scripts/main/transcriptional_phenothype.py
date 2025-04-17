# import muon.pp.filter_obs as filter_obs ???
import muon as mu
import scanpy as sc
import argparse
import pandas as pd
import sys
from anndata import AnnData

sys.path.append("/opt/scripts/utility")
from helpers import update_validation_metrics



def compute_mmc_taxonomy(mmc_taxonomy_file: str) -> pd.DataFrame:
    """
    Load the MMC taxonomy file
    """
    return pd.read_csv(mmc_taxonomy_file)

def find_high_fidelity_mappings(mmc_taxonomy_file: str) -> pd.DataFrame:
    """
    Find the high-fidelity mappings in the MMC taxonomy file
    """
    return pd.read_csv(mmc_taxonomy_file)




def main(args: argparse.Namespace):
    """
    basic logic with args as input

    """


    # MMC
    mmc_taxonomy = compute_mmc_taxonomy(args.mmc_taxonomy_file)


    # filter types for "class" by correlation and bootstrap probability > 0.5
    # assign the rest of the cells to "unknown"
    high_fidelity_mappings = find_high_fidelity_mappings(mmc_taxonomy)

    # save the filtered adata
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
