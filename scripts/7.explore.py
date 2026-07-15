"""
Does not respect overwrite.

Loads the per-gene feature-selected profiles output by 6.feature-select
(gene-level only, not guide-level) and plots a UMAP embedding where each
point is a gene, with the negative control gene(s) highlighted.

Outputs the following files:
- data_exploration/{plate}_{data_set_name}_gene_normalized_{group}_UMAP.png
- data_exploration/{plate}_{data_set_name}_gene_normalized_{group}_feature_selected_mAP.csv.gz
- data_exploration/{plate}_{data_set_name}_gene_normalized_{group}_feature_selected_mAP_scatterplot.png
- data_exploration/{plate}_{data_set_name}_gene_normalized_{group}_feature_selected_{gene}_cosine_similarity.csv.gz
"""

import os
import sys
import glob
import logging
import traceback
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import umap
from copairs import map
from sklearn.metrics.pairwise import cosine_similarity

recipe_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(recipe_path, "utils"))
from io_utils import load_configs, read_csvs_with_chunksize, parse_data_set

# Configure logging
logfolder = os.path.join(recipe_path, "logs")
if not os.path.isdir(logfolder):
    os.mkdir(logfolder)
logging.basicConfig(
    filename=os.path.join(logfolder, "7.explore.log"),
    level=logging.INFO,
)


def printandlog(msg, type="info"):
    print(msg)
    if type == "warning":
        logging.warning(msg)
    else:
        logging.info(msg)


def handle_excepthook(exc_type, exc_value, exc_traceback):
    logging.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))
    traceback_details = "\n".join(traceback.extract_tb(exc_traceback).format())
    print(f"Uncaught Exception: {traceback_details}")


sys.excepthook = handle_excepthook


def make_heatmap(groupdf, cols, outpath):
    corr = groupdf.set_index(cols)[
        [x for x in groupdf.columns if "Metadata" not in x]
    ].T.corr()
    # Set up the matplotlib figure
    fig, ax = plt.subplots(figsize=(11, 9))

    # Generate a custom diverging colormap
    sns.heatmap(
        corr,
        cmap="vlag",
        vmax=0.3,
        center=0,
        square=True,
        linewidths=0.5,
        cbar_kws={"shrink": 0.5},
        ax=ax,
    )
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)


def explore(path_to_defaults_config, path_to_experiment_config):
    printandlog("Starting 7.explore.")
    defaults_config, experiment_config = load_configs(
        path_to_defaults_config, path_to_experiment_config
    )

    # define experiment variables
    control_genes = experiment_config["control_genes"]
    poscon_genes = experiment_config["poscon_genes"]
    data_sets = experiment_config["data_sets"]

    # define defaults variables
    out_root = defaults_config["directory_structure"]["root"]
    profiles_folder = defaults_config["directory_structure"]["profiles"]
    exploration_folder = defaults_config["directory_structure"]["exploration"]
    gene_col = defaults_config["process"]["process_SBS"]["gene_col"]
    guide_col = defaults_config["process"]["process_SBS"]["barcode_col"]
    group = defaults_config["process"]["feature_select"]["group"]

    profiles_out = os.path.join(out_root, profiles_folder)
    outdir_exploration = os.path.join(out_root, exploration_folder)
    if not os.path.exists(outdir_exploration):
        os.makedirs(outdir_exploration, exist_ok=True)

    for data_set_name in data_sets.keys():
        printandlog(f"Starting {data_set_name}")
        batch, plate_list, list_plate_well_tuples = parse_data_set(
            data_sets[data_set_name]
        )
        for plate in plate_list:
            for agg in ["gene", "guide"]:
                groupdf = read_csvs_with_chunksize(
                    os.path.join(
                        profiles_out,
                        f"{plate}_{data_set_name}_{agg}_normalized_feature_selected_by{group}.csv.gz",
                    )
                )
                #
                # Create UMAP embedding - gene and guide maps
                #

                # create color dictionary for consistency in plot
                all_unique_genes = groupdf[f"Metadata_{gene_col}"].unique()
                colors = sns.color_palette("husl", len(all_unique_genes))
                gene_palette = dict(zip(all_unique_genes, colors))

                feature_cols = [x for x in groupdf.columns if "Metadata" not in x]

                embedding = umap.UMAP().fit_transform(groupdf[feature_cols])
                umap_df = pd.DataFrame(embedding, columns=["UMAP1", "UMAP2"])
                umap_df[f"Metadata_{gene_col}"] = groupdf[f"Metadata_{gene_col}"].values

                fig, ax = plt.subplots()
                sns.scatterplot(
                    data=umap_df,
                    x="UMAP1",
                    y="UMAP2",
                    hue=f"Metadata_{gene_col}",
                    palette=gene_palette,
                    ax=ax,
                )
                sns.scatterplot(
                    data=umap_df.loc[
                        umap_df[f"Metadata_{gene_col}"].isin(control_genes)
                    ],
                    x="UMAP1",
                    y="UMAP2",
                    marker="h",
                    edgecolor="black",
                    hue=f"Metadata_{gene_col}",
                    palette=gene_palette,
                    legend=False,
                    ax=ax,
                )
                if poscon_genes:
                    for gene in poscon_genes:
                        sns.scatterplot(
                            data=umap_df.loc[umap_df[f"Metadata_{gene_col}"] == gene],
                            x="UMAP1",
                            y="UMAP2",
                            edgecolor="black",
                            hue=f"Metadata_{gene_col}",
                            palette=gene_palette,
                            marker="*",
                            legend=False,
                            ax=ax,
                        )
                ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1))
                ax.set_title(f"{plate} {data_set_name} {agg} normalized by {group}")
                fig.savefig(
                    os.path.join(
                        outdir_exploration,
                        f"{plate}_{data_set_name}_{agg}_normalized_{group}_UMAP.png",
                    ),
                    dpi=300,
                    bbox_inches="tight",
                )
                plt.close(fig)

                if poscon_genes:
                    for gene in poscon_genes:
                        #
                        # Calculate cosine similarity to positive controls - gene and guide
                        #
                        cov = cosine_similarity(
                            groupdf[
                                [x for x in groupdf.columns if "Metadata" not in x]
                            ].values,
                            groupdf.loc[groupdf[f"Metadata_{gene_col}"] == gene, :][
                                [x for x in groupdf.columns if "Metadata" not in x]
                            ].values,
                        )
                        cos_sim_df = pd.DataFrame(
                            cov,
                            columns=["cosine_similarity"],
                            index=groupdf[f"Metadata_{gene_col}"],
                        )
                        cos_sim_df.to_csv(
                            os.path.join(
                                outdir_exploration,
                                f"{plate}_{data_set_name}_{agg}_normalized_{group}_feature_selected_{gene}_cosine_similarity.csv.gz",
                            ),
                            index=True,
                        )

                #
                # Create heatmap - gene and guide maps
                #
                if agg == "gene":
                    make_heatmap(
                        groupdf,
                        [f"Metadata_{gene_col}"],
                        os.path.join(
                            outdir_exploration,
                            f"{plate}_{data_set_name}_gene_normalized_{group}_UMAP.png",
                        ),
                    )
                if agg == "guide":
                    make_heatmap(
                        groupdf,
                        [f"Metadata_{gene_col}", f"Metadata_{guide_col}"],
                        os.path.join(
                            outdir_exploration,
                            f"{plate}_{data_set_name}_guide_normalized_{group}_UMAP.png",
                        ),
                    )

            #
            # Calculate mAP for phenotypic activity - guide only
            #
            df_activity = groupdf.copy(deep=True).dropna()
            df_activity = df_activity.loc[
                :, df_activity.nunique() > 1
            ]  # remove constant columns

            # make default value equal to row index
            df_activity["Metadata_reference_index"] = df_activity.index
            # make index equal to -1 for all treatment replicates (non-DMSO)
            df_activity.loc[
                ~(df_activity[f"Metadata_{gene_col}"].isin(control_genes)),
                "Metadata_reference_index",
            ] = -1
            df_activity.insert(
                0,
                "Metadata_reference_index",
                df_activity.pop("Metadata_reference_index"),
            )

            # positive pairs are replicates of the same treatment
            # negative pairs are replicates of different treatments
            pos_sameby = [f"Metadata_{gene_col}", "Metadata_reference_index"]
            pos_diffby = []
            neg_sameby = []
            neg_diffby = [f"Metadata_{gene_col}", "Metadata_reference_index"]

            metadata = df_activity.filter(regex="^Metadata")
            profiles = df_activity.filter(regex="^(?!Metadata)").values

            replicate_aps = map.average_precision(
                metadata, profiles, pos_sameby, pos_diffby, neg_sameby, neg_diffby
            )
            replicate_maps = map.mean_average_precision(
                replicate_aps, pos_sameby, null_size=10000, threshold=0.05, seed=0
            )
            replicate_maps["-log10(p-value)"] = -replicate_maps[
                "corrected_p_value"
            ].apply(np.log10)

            # Output mAP results, sorted by mAP
            df = replicate_maps.reset_index().sort_values(by="mean_average_precision")
            df.to_csv(
                os.path.join(
                    outdir_exploration,
                    f"{plate}_{data_set_name}_guide_normalized_{group}_feature_selected_mAP.csv.gz",
                ),
                index=False,
            )

            # Plot mAP results
            fig, ax = plt.subplots()
            sns.scatterplot(
                data=df,
                x="mean_average_precision",
                y="-log10(p-value)",
                ax=ax,
            )
            if poscon_genes:
                legend_elements = []
                for gene in poscon_genes:
                    sns.scatterplot(
                        data=df.loc[df[f"Metadata_{gene_col}"] == gene],
                        x="mean_average_precision",
                        y="-log10(p-value)",
                        edgecolor="black",
                        marker="^",
                        legend=False,
                        ax=ax,
                    )
                    legend_elements += [
                        Line2D(
                            [0],
                            [0],
                            marker="^",
                            color="tab:blue",
                            label=gene,
                            markersize=5,
                            lw=0,
                        )
                    ]
                ax.legend(handles=legend_elements, loc="best")
            ax.set_title(f"{plate} {data_set_name} guide normalized by {group}")
            outpath = os.path.join(
                outdir_exploration,
                f"{plate}_{data_set_name}_guide_normalized_{group}_feature_selected_mAP_scatterplot.png",
            )
            fig.savefig(outpath, dpi=300, bbox_inches="tight")
            plt.close(fig)

    printandlog("Finished 7.explore.")


#################################
# MAIN USER INTERACTION
#################################

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Pass paths to defaults config and experiment config.")
        sys.exit()
    path_to_defaults_config = sys.argv[1]
    path_to_experiment_config = sys.argv[2]
    explore(path_to_defaults_config, path_to_experiment_config)
