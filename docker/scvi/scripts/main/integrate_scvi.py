import os
import argparse
import scvi
import anndata as ad

def integrate_with_scvi(adata: ad.AnnData, batch_key: str, latent_key: str) -> tuple[ad.AnnData, scvi.model.SCVI]:
    """
    fit scvi model to AnnData object
    """
    
    
    ## fixed parameters
    n_latent = 30
    n_layers = 1
    train_size = 0.85
    scvi_epochs = 200
    batch_size = 1024
    # accelerator = 'gpu'
    # devide = "cuda:0"
    dispersion = "gene"
    gene_likelihood = "zinb"
    latent_distribution = "normal"
    early_stopping = True
    ## DEPRICATE these training parameters. defaults are good
    # plan_kwargs = {"lr_factor": 0.1, "lr_patience": 20, "reduce_lr_on_plateau": True}
    # early_stopping_patience = 25

  
    # integrate the data with `scVI`
    noise = ["doublet_score", "pct_counts_mt", "pct_counts_rb"]
    categorical_covariate_keys = None
    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        batch_key=batch_key,
        continuous_covariate_keys=noise,
        categorical_covariate_keys=categorical_covariate_keys,
    )

    model = scvi.model.SCVI(
        adata,
        n_layers=n_layers,
        n_latent=n_latent,
        dispersion=dispersion,
        gene_likelihood=gene_likelihood,
    )

    model.train(
        train_size=train_size,
        max_epochs=scvi_epochs,
        early_stopping=early_stopping,
        # early_stopping_patience=early_stopping_patience,
        # plan_kwargs=plan_kwargs,
    )

    adata.obsm[latent_key] = model.get_latent_representation()  # type: ignore
                        
    return adata, model


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run scVI integration")
    parser.add_argument(
        "--latent-key",
        dest="latent_key",
        type=str,
        default="X_scvi",
        help="Latent key to save the scvi latent to",
    )
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
        "--output-scvi-dir",
        dest="output_scvi_dir",
        type=str,
        help="Output folder to save `scvi` model",
    )

    # TODO: optional scvi arguments
    args = parser.parse_args()

    # 0. load data
    adata = ad.read_h5ad(args.adata_input)  # type: ignore

    # 2. process data
    adata,model = integrate_with_scvi(adata, args.batch_key, args.latent_key)
    

    # 3. save the integrated adata and scvi model
    model.save(args.output_scvi_dir, overwrite=True)
    adata.write_h5ad(filename=args.adata_output, compression="gzip")
