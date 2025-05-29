#!/bin/bash

ROOT="/home/jupyter/workspace/ws_files/data/preprocess/"
# Set variables
COHORT_ID="pmdbs_sc_rnaseq_cohort_analysis_team-scherzer"
BATCH_KEY="sample"
LABEL_KEY="cell_type"
N_TOP_GENES=5000

# Input/output paths
INPUT_PATH="${ROOT}/data/h5_files"  # Path containing .h5 files
PREPROCESSED_ADATA_DIR="${ROOT}/data/preprocessed_adata"

# Filter and normalize
python filter.py \
    --adata-input "${PREPROCESSED_ADATA_DIR}/${COHORT_ID}.merged_adata_object.h5ad" \
    --adata-output "${PREPROCESSED_ADATA_DIR}/${COHORT_ID}.filtered.h5ad" \
    --output-validation-file "${PREPROCESSED_ADATA_DIR}/${COHORT_ID}.validation_metrics.csv"

# now do mmc
python mmc.py \
    --adata-input "${PREPROCESSED_ADATA_DIR}/${COHORT_ID}.filtered.h5ad" \
    --output-name "${PREPROCESSED_ADATA_DIR}/${COHORT_ID}" \
    --mmc-taxonomy-path "${PREPROCESSED_ADATA_DIR}/precomputed_stats.20231120.sea_ad.MTG.h5" 
    # \
    # --mmc-out-path "$HOME/Projects/ASAP/pmdbs-sc-rnaseq-wf"

python process.py \
    --adata-input "${PREPROCESSED_ADATA_DIR}/${COHORT_ID}.filtered.h5ad" \
    --batch-key "${BATCH_KEY}" \
    --adata-output "${COHORT_ID}.filtered_normalized.h5ad" \
    --n-top-genes "${N_TOP_GENES}" \
    --output-all-genes "${PREPROCESSED_ADATA_DIR}/${COHORT_ID}.all_genes.csv" \
    --output-hvg-genes "${PREPROCESSED_ADATA_DIR}/${COHORT_ID}.hvg_genes.csv" \
    --output-validation-file "${PREPROCESSED_ADATA_DIR}/${COHORT_ID}.final_validation_metrics.csv"

    # --marker-genes "${CELL_TYPE_MARKERS}" \

# Transcriptional phenotype
python  transcriptional_phenotype.py \
    --adata-input "${PREPROCESSED_ADATA_DIR}/${COHORT_ID}.filtered_normalized.h5ad" \
    --mmc-results "${PREPROCESSED_ADATA_DIR}/${COHORT_ID}.mmc.SEAAD_results.csv" \
    --output-cell-types-file "${PREPROCESSED_ADATA_DIR}/${COHORT_ID}.mmc.cell_types.parquet" \
    --adata-output "${PREPROCESSED_ADATA_DIR}/${COHORT_ID}.annotated.h5ad"
#"${COHORT_ID}.mmc_results.csv" \
    
# Integration with scVI
python integrate_scvi.py \
    --adata-input "${PREPROCESSED_ADATA_DIR}/${COHORT_ID}.annotated.h5ad" \
    --batch-key "${BATCH_KEY}" \
    --adata-output "${PREPROCESSED_ADATA_DIR}/${COHORT_ID}.integrated.h5ad" \
    --output-scvi-dir "${PREPROCESSED_ADATA_DIR}/${COHORT_ID}.scvi_model" \
    --output-scanvi-dir "${PREPROCESSED_ADATA_DIR}/${COHORT_ID}.scanvi_model" \
    --output-cell-types-file "${PREPROCESSED_ADATA_DIR}/${COHORT_ID}.scanvi.cell_types.parquet"











# Clustering and UMAP
python3 /opt/scripts/main/clustering_umap.py \
    --adata-input "${COHORT_ID}.integrated.h5ad" \
    --adata-output "${COHORT_ID}.clustered.h5ad"

# Add Harmony integration
python3 /opt/scripts/main/add_harmony.py \
    --adata-input "${COHORT_ID}.clustered.h5ad" \
    --batch-key "${BATCH_KEY}" \
    --adata-output "${COHORT_ID}.final.h5ad" \
    --output-metadata-file "${COHORT_ID}.final_metadata.csv"

# Calculate artifact metrics
python3 /opt/scripts/main/artifact_metrics.py \
    --adata-input "${COHORT_ID}.final.h5ad" \
    --batch-key "${BATCH_KEY}" \
    --label-key "${LABEL_KEY}" \
    --output-report-dir "scib_report"

# Plot features and groups
python3 /opt/scripts/main/plot_feats_and_groups.py \
    --adata-input "${COHORT_ID}.final.h5ad"
