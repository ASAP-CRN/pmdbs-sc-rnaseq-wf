#!/bin/bash

ROOT="/home/ergonyc/Projects/ASAP/pmdbs-sc-rnaseq-wf"

# Set variables
COHORT_ID="pmdbs_sc_rnaseq_cohort_analysis_team-scherzer"
BATCH_KEY="sample"
LABEL_KEY="cell_type"
N_TOP_GENES=5000

# Input/output paths
ROOT="/home/ergonyc/Projects/ASAP/pmdbs-sc-rnaseq-wf"
INPUT_PATH="${ROOT}/data/h5_files"  # Path containing .h5 files
PREPROCESSED_ADATA_DIR="${ROOT}/data/preprocessed_adata"

# Create output directory if it doesn't exist
mkdir -p ${PREPROCESSED_ADATA_DIR}

# Apply prep_metadata.py to all .h5 files in INPUT_PATH
for h5_file in ${INPUT_PATH}/*.h5; do
    filename=$(basename -- "$h5_file")
    sample_id="${filename%.*}"  # Remove extension to get sample_id
    
    echo "Processing $h5_file..."
    
    /home/ergonyc/mambaforge/envs/scib-nb-cuda12.4/bin/python /opt/scripts/main/prep_metadata.py \
        --adata-input "$h5_file" \
        --sample-id "$sample_id" \
        --batch "batch1" \
        --dataset "$COHORT_ID" \
        --team "team-scherzer" \
        --adata-output "${PREPROCESSED_ADATA_DIR}/${sample_id}.preprocessed.h5ad"
done

# Create a file listing all preprocessed adata paths for merge_and_plot_qc.py
find ${PREPROCESSED_ADATA_DIR} -name "*.preprocessed.h5ad" | awk '{print $0"\t"$0}' > "${PREPROCESSED_ADATA_DIR}/adata_samples_paths.tsv"

# Merge and plot QC metrics
python3 /opt/scripts/main/merge_and_plot_qc.py \
    --adata-objects-fofn "${PREPROCESSED_ADATA_DIR}/adata_samples_paths.tsv" \
    --adata-output "${COHORT_ID}.merged_adata_object.h5ad" \
    --output-metadata-file "${COHORT_ID}.initial_metadata.csv" \
    --output-validation-file "${COHORT_ID}.validation_metrics.csv"

# Filter and normalize
/home/ergonyc/mambaforge/envs/scib-nb-cuda12.4/bin/python  /opt/scripts/main/filter.py \
    --adata-input "${COHORT_ID}.merged_adata_object.h5ad" \
    --adata-output "${COHORT_ID}.filtered.h5ad" \
    --output-validation-file "${COHORT_ID}.validation_metrics.csv"

# now do mmc
/home/ergonyc/mambaforge/envs/mmc_base/bin/python  /opt/scripts/main/mmc.py \
    --adata-input "${COHORT_ID}.filtered.h5ad" \
    --output-name "${COHORT_ID}" \
    --mmc-taxonomy-path "precomputed_stats.20231120.sea_ad.MTG.h5" 
    # \
    # --mmc-out-path "$HOME/Projects/ASAP/pmdbs-sc-rnaseq-wf"


/home/ergonyc/mambaforge/envs/scib-nb-cuda12.4/bin/python  /opt/scripts/main/process.py \
    --adata-input "${COHORT_ID}.filtered.h5ad" \
    --batch-key "${BATCH_KEY}" \
    --adata-output "${COHORT_ID}.filtered_normalized.h5ad" \
    --n-top-genes "${N_TOP_GENES}" \
    --output-all-genes "${COHORT_ID}.all_genes.csv" \
    --output-hvg-genes "${COHORT_ID}.hvg_genes.csv" \
    --output-validation-file "${COHORT_ID}.final_validation_metrics.csv"

    # --marker-genes "${CELL_TYPE_MARKERS}" \

# Transcriptional phenotype
/home/ergonyc/mambaforge/envs/scib-nb-cuda12.4/bin/python  /opt/scripts/main/transcriptional_phenotype.py \
    --adata-input "${COHORT_ID}.filtered_normalized.h5ad" \
    --mmc-results "/home/ergonyc/tmp/pmdbs_sc_rnaseq_cohort_analysis_team-scherzer.mmc.SEAAD_results.csv" \
    --output-cell-types-file "${COHORT_ID}.mmc.cell_types.parquet" \
    --adata-output "${COHORT_ID}.annotated.h5ad"
#"${COHORT_ID}.mmc_results.csv" \
    
# Integration with scVI
/home/ergonyc/mambaforge/envs/scib-nb-cuda12.4/bin/python  /opt/scripts/main/integrate_scvi.py \
    --adata-input "${COHORT_ID}.annotated.h5ad" \
    --batch-key "${BATCH_KEY}" \
    --adata-output "${COHORT_ID}.integrated.h5ad" \
    --output-scvi-dir "${COHORT_ID}.scvi_model" \
    --output-scanvi-dir "${COHORT_ID}.scanvi_model" \
    --output-cell-types-file "${COHORT_ID}.scanvi.cell_types.parquet"











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
