# pip install -U git+https://github.com/AllenInstitute/cell_type_mapper.git

from cell_type_mapper.cli.map_to_on_the_fly_markers import OnTheFlyMapper
from pathlib import Path
import os
import argparse
# works best on CPU
os.environ["AIBS_BKP_USE_TORCH"] = "false"

# os.environ["NUMEXPR_NUM_THREADS"] = "1"
# os.environ["MKL_NUM_THREADS"] = "1"
# os.environ["OMP_NUM_THREADS"] = "1"

HOME = Path.home()
TMP_DIR = HOME / "tmp"
if not TMP_DIR.exists():
    TMP_DIR.mkdir()

# these wer some initial defaults
CHUNK_SIZE = 40000
N_RUNNERS_UP = 5
RNG_SEED = 11235813
N_PROCESSORS = 8
MAX_GB = 48.0


##  how do we want to deal with this resource... it just needs to be downloaded once and then cached
##  the probably just copy it to the local directory first.
## TODO: make sure https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/mapmycells/SEAAD/20240831/precomputed_stats.20231120.sea_ad.MTG.h5
# can also use this: !pip install -U git+https://github.com/alleninstitute/abc_atlas_access >& data/scratch/junk.txt

PRECOMPUTED_STATS = (
    f"examples/data/abc_atlas_data/precomputed_stats.20231120.sea_ad.MTG.h5"
)



def main(args: argparse.Namespace):
    """
    basic logic with args as input

    """

    # not sure if we should use a tmp dir or not.
    tmp_dir = TMP_DIR 
    if not tmp_dir.exists():
        tmp_dir.mkdir()

    REFERENCE = "SEAAD"

    file_root = args.result_name

    EXTENDED_RESULTS = f"{file_root}.mmc.{REFERENCE}.json"
    LOG_FILE = f"{file_root}.mmc.{REFERENCE}_log.txt"
    CSV_RESULTS = f"{file_root}.mmc.{REFERENCE}_results.csv"


    # MMC
    config = {
        # "query_path": f"{ASAP_DATA}/{FILE_ROOT}.h5ad",
        "query_path": args.adata_input,
        "tmp_dir": f"{TMP_DIR}",
        "extended_result_path": f"{tmp_dir / EXTENDED_RESULTS}",
        "csv_result_path": f"{tmp_dir / CSV_RESULTS}",
        "log_path": f"{tmp_dir / LOG_FILE}",
        "cloud_safe": False,
        "verbose_csv": True,
        "n_processors": N_PROCESSORS,
        "max_gb": MAX_GB,
        "map_to_ensembl": True,
        "type_assignment": {
            "normalization": "raw",
            "bootstrap_iteration": 100,
            "bootstrap_factor": 0.5,
            "chunk_size": CHUNK_SIZE,
            "n_runners_up": N_RUNNERS_UP,
            "rng_seed": RNG_SEED,
        },
        "precomputed_stats": {"path": args.mmc_taxonomy_path},
        "reference_markers": {"log2_fold_min_th": 0.5},
        "query_markers": {"n_per_utility": 15, "genes_at_a_time": 1},
    }

    runner = OnTheFlyMapper(args=[], input_data=config)
    runner.run()

    # this creates the results files. .json, .csv, .txt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Transcriptional Phenotype (MMC)")
    parser.add_argument(
        "--adata-input",
        dest="adata_input",
        type=str,
        help="AnnData object for a dataset",
    )

    parser.add_argument(
        "--mmc-taxonomy-path",
        dest="mmc_taxonomy_path",
        type=str,
        help="Path to MMC taxonomy files",
        default=f"{PRECOMPUTED_STATS}"
    )
 
    parser.add_argument(
        "--output-name",
        dest="result_name",
        type=str,
        help="Output files to write results to",
    )

    args = parser.parse_args()
    main(args)


# # SEA-AD "references"
# MTG
# DLPFC
# MTG_DLPFC

# look at discrepancies.

# ## xylena reference...
# xylena_train for Taxonomy

# test on xylena_test + xylena_query
