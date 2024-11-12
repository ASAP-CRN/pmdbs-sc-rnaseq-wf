import argparse
import scanpy as sc


def make_f_and_g_plots(adata: sc.AnnData, group_name: str, feature_name: str):
    # Set CPUs to use for parallel computing
    sc._settings.ScanpyConfig.n_jobs = -1

    # Set working directory and load packages
    sc.settings.verbosity = 1
    sc.settings.figdir = "plots/"
    sc.settings.set_figure_params(
        dpi=100, fontsize=10, dpi_save=300, format="png", figsize=("12", "8")
    )  # type: ignore

    features_list = args.feature.split(",")
    plot_features = [x for x in features_list if x in adata.obs.columns]
    file_name = group_name + "_features_umap.png"
    sc.pl.embedding(
        adata,
        basis="umap",
        color=plot_features,
        frameon=False,
        show=False,
        ncols=1,
        save=file_name,
    )

    groups_list = args.group.split(",")
    plot_groups = [x for x in groups_list if x in adata.obs.columns]
    file_name = feature_name + "_groups_umap.png"
    sc.pl.embedding(
        adata,
        basis="umap",
        color=plot_groups,
        frameon=False,
        show=False,
        ncols=1,
        save=file_name,
    )


def main(args: argparse.Namespace):
    """
    basic logic with args as input

    """
    adata = sc.read_h5ad(args.adata_input, backed="r")  # type: ignore
    group_name = args.output_group_umap_plot_prefix
    feature_name = args.output_feature_umap_plot_prefix
    make_f_and_g_plots(adata, group_name, feature_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot groups")
    parser.add_argument(
        "--adata-input",
        dest="adata_input",
        type=str,
        help="AnnData object for a dataset",
    )
    parser.add_argument(
        "--group", dest="group", type=str, help="Group to plot umaps for"
    )
    parser.add_argument(
        "--output-group-umap-plot-prefix",
        dest="output_group_umap_plot_prefix",
        type=str,
        help="Output file prefix to write the group umap plot to; will write both a PDF and a PNG.",
    )
    parser.add_argument(
        "--feature", dest="feature", type=str, help="Feature to plot umaps for"
    )
    parser.add_argument(
        "--output-feature-umap-plot-prefix",
        dest="output_feature_umap_plot_prefix",
        type=str,
        help="Output file to write the feature umap plot to",
    )

    args = parser.parse_args()
    main(args)
