# scVI analysis Pipeline python scripts


### _PREPROCESS_
- _pre-preprocessing_: executed by wdl [`cellbender :: remove_technical_artifacts`](../../../workflows/preprocess/preprocess.wdl#L327-L350)
- _doublets_ + _metadata_: [`prep_metadata.py`](./main/prep_metadata.py), calculates `scrublet` metrics and adds additional metadata.


### _MERGE_AND_QC: merge and plot QC metrics_
- [`merge_and_plot_qc.py`](./main/merge_and_plot_qc.py)
    - Merges adatas.
    - plot general QC metrics for all cells
    - Dump final metadata

    - NOTE:  merging here could be inefficient.  Stubs of code to create gene_list to pre-filter.
        - Consider filtering cells first... 'processing' filters genes to "highly variable".


### _QC: FILTER_ 
- _filtering_: [`filter.py`](./main/filter.py)


### _FEATURE SELECTION_ 
- _processing_: [`process.py`](./main/process.py)
    - Normalize + feature selection (i.e. identification of highly variable genes).
    - add PCA

### _INTEGRATION_
- [`integrate_scvi.py`](./main/integrate_scvi.py)
    - `scVI` integration.


### _CLUSTERING_
- _clustering_: `umap` ([`clustering_umap.py`](./main/clustering_umap.py))
    - updated to do leiden at 4 resolutions - [0.05, 0.1, 0.2, 0.4]
    - In the future we may choose `mde` (`clustering_mde.py`) over `umap`, as it is super fast and efficient on a GPU, and the embeddings are only useful for visualization so the choice is semi-arbitrary.


### _ANNOTATION_
- NOTE:  considering that this might be an "additional curation" item for cross-team harmonized artifacts.  This annotation should probably be done separately on different "atlas" level subsets of the data.
- [`annotate_cells.py`](./main/annotate_cells.py).  Use cellassign and a list of marker genes.  Need to curate a list of marker genes apropriate for all the brain regions we are processing.  Currently using CARD cortical list of genes.  NOTE: this is not annotating the "clusters" but the cells based on marker gene expression.


### _HARMONY INTEGRATION_
- [`add_harmony.py`](./main/add_harmony.py)
    - Add and Harmony integration.
    - Dump final metadata


### _`SCIB` METRICS_
- [`artifact_metrics.py`](./main/artifact_metrics.py)
    - Compute `scib` metrics on final artifacts and generate a report.


### _PLOTTING_
- [`plot_feats_and_groups.py`](./main/plot_feats_and_groups.py).  
- features: 'n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_rb', 'doublet_score','S_score','G2M_score'
- groups: "sample", "batch", "cell_type"


### CONSIDERATIONS.
The "best" tools are probably `scAR`+`SOLO`+`scVI`.  Tradeoffs are probably some computation.  However, without further testing we will default to `cellbender`+`scrublet` for removing ambient and creating the doublet metrics. `scVI` to be used for integration.

