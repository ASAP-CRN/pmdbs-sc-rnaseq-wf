#!/bin/bash

ROOT="/home/jupyter/workspace/ws_files/data/preprocess"
# Set variables
COHORT_ID="pmdbs_sc_rnaseq_cohort_analysis_team-scherzer"
BATCH_KEY="sample"
LABEL_KEY="cell_type"
N_TOP_GENES=5000

# Input/output paths
INPUT_PATH="${ROOT}"  # Path containing .h5 files
PREPROCESSED_ADATA_DIR="${ROOT}/data/preprocessed_adata"

# Create output directory if it doesn't exist
mkdir -p ${PREPROCESSED_ADATA_DIR}

# # Apply prep_metadata.py to all .h5 files in INPUT_PATH
# for h5_file in ${INPUT_PATH}/*.h5; do
#     filename=$(basename -- "$h5_file")
#     sample_id="${filename%.*}"  # Remove extension to get sample_id
    
#     echo "Processing $h5_file..."
    
#     python3 prep_metadata.py \
#         --adata-input "$h5_file" \
#         --sample-id "$sample_id" \
#         --batch "batch1" \
#         --dataset "$COHORT_ID" \
#         --team "team-scherzer" \
#         --adata-output "${PREPROCESSED_ADATA_DIR}/${sample_id}.preprocessed.h5ad"
# done

# Create a file listing all preprocessed adata paths for merge_and_plot_qc.py
find ${INPUT_PATH} -name "*.preprocessed.h5ad" | awk '{print $0"\t"$0}' > "${PREPROCESSED_ADATA_DIR}/adata_samples_paths.tsv"

# Merge and plot QC metrics
python3 merge_and_plot_qc.py \
    --adata-objects-fofn "${PREPROCESSED_ADATA_DIR}/adata_samples_paths.tsv" \
    --adata-output "${PREPROCESSED_ADATA_DIR}/${COHORT_ID}.merged_adata_object.h5ad" \
    --output-metadata-file "${PREPROCESSED_ADATA_DIR}/${COHORT_ID}.initial_metadata.csv" \
        --output-validation-file "${PREPROCESSED_ADATA_DIR}/${COHORT_ID}.validation_metrics.csv"
