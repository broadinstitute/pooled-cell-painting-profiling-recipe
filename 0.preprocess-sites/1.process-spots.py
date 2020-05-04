"""
1.process-spots.py

To view and change parameters used in this file, see `process-spots` in the config.

We capture several key observations for each site in a pooled cell painting experiment.
Each site is a unique image taken by the microscope.

The observations are:

1. Key statistics
    * Number of spots
    * Number of cells
    * Number of unique genes
2. Summary figures
    * Distribution of barcode scores
    * Distribution of number of barcodes per cell
    * Relationship of number of barcodes to average score
    * Several figures for gene scores
        * Distribution of number of high quality cells per gene
3. Results tables
    * Cells with quality categories, gene infection, mean scores, and barcode counts
    * Cells with quality categories, CRISPR infection by guide, mean scores and counts
    * Summary of cell count qualities
    * Number of different quality cells per gene

We assign cell quality estimates to each cell.

All sites are processed and the results are saved in site specific folders
"""

import os
import sys
import pathlib
import argparse
import pandas as pd
import yaml

import plotnine as gg
import matplotlib.pyplot as plt
import seaborn as sns

from scripts.spot_utils import (
    spot_counts_per_cell_histogram,
    spot_score_histogram,
    spot_count_score_jointplot,
    category_counts,
)

sys.path.append(os.path.join("..", "scripts"))
from config_utils import get_config
from cell_quality_utils import CellQuality

parser = argparse.ArgumentParser()
parser.add_argument(
    "--config_file",
    help="configuration yaml file for preprocessing pipeline",
    default="site_processing_config.yaml",
)
args = parser.parse_args()
config_file = args.config_file

config = get_config(config_file)

core_args = config["core"]
spot_args = config["process-spots"]

project = core_args["project"]
batch = core_args["batch"]
batch_dir = core_args["batch_dir"]
quality_func = core_args["cell_quality"]
control_barcodes = core_args["control_barcodes"]
trash_files = core_args["trash_files"]

id_cols = core_args["id_cols"]
spot_parent_cols = core_args["parent_cols"]["spots"]

output_basedir = spot_args["output_basedir"]
barcode_cols = spot_args["barcode_cols"]
gene_cols = spot_args["gene_cols"]
location_cols = spot_args["location_cols"]
spot_score_cols = spot_args["spot_score_cols"]
foci_cols = spot_args["foci_cols"]
cell_filter = spot_args["cell_filter"]

barcode_foci_cols = id_cols + location_cols + spot_parent_cols
all_foci_cols = list(
    set(
        id_cols + location_cols + foci_cols + spot_score_cols + barcode_cols + gene_cols
    )
)

cell_quality = CellQuality(quality_func)
cell_category_dict = cell_quality.define_cell_quality()
cell_category_df = pd.DataFrame(cell_category_dict, index=["Cell_Class"])

sites = [x for x in os.listdir(batch_dir) if x not in trash_files]
num_sites = len(sites)

for site in sites:
    print(f"Now processing {site}...")
    output_dir = pathlib.PurePath(f"{output_basedir}/{batch}/spots/{site}")
    os.makedirs(output_dir, exist_ok=True)

    # Load spot data
    try:
        barcode_file = pathlib.PurePath(f"{batch_dir}/{site}/BarcodeFoci.csv")
        barcodefoci_df = pd.read_csv(barcode_file)

        foci_file = pathlib.PurePath(f"{batch_dir}/{site}/Foci.csv")
        foci_df = pd.read_csv(foci_file)
    except FileNotFoundError:
        print(f"{site} data not found")

    image_number = foci_df.ImageNumber.unique()[0]

    try:
        # Confirm that image number and object number are aligned
        pd.testing.assert_frame_equal(
            barcodefoci_df.loc[:, id_cols], foci_df.loc[:, id_cols], check_names=True
        )

        pd.testing.assert_frame_equal(
            barcodefoci_df.loc[:, location_cols],
            foci_df.loc[:, location_cols],
            check_names=True,
        )
    except AssertionError:
        print(f"{site} data not aligned between foci files")

    # Merge spot data files
    complete_foci_df = barcodefoci_df.loc[:, barcode_foci_cols].merge(
        foci_df.loc[:, all_foci_cols],
        left_on=id_cols + location_cols,
        right_on=id_cols + location_cols,
        how="inner",
    )

    null_spot_df = complete_foci_df.query("Parent_Cells == 0")
    cell_spot_df = complete_foci_df.query("Parent_Cells != 0")

    num_assigned_cells = len(cell_spot_df.Parent_Cells.unique())
    num_unassigned_spots = null_spot_df.shape[0]
    num_assigned_spots = cell_spot_df.shape[0]

    # Figure 1 - histogram of barcode counts per cell
    fig_file = pathlib.PurePath(f"{output_dir}/num_spots_per_cell_histogram.png")
    spot_counts_per_cell_histogram(cell_spot_df, spot_parent_cols, fig_file)

    # Figure 2 - histogram of barcode scores per spot
    fig_file = pathlib.PurePath(f"{output_dir}/barcode_scores_per_spot_histogram.png")
    spot_score_histogram(cell_spot_df, spot_score_cols, fig_file)

    # Figure 3 - Joint plot of relationship of barcode counts per cell and mean score
    fig_file = pathlib.PurePath(
        f"{output_dir}/per_cell_barcode_count_by_mean_score_jointplot.png"
    )
    spot_count_score_jointplot(
        cell_spot_df, spot_parent_cols[0], spot_score_cols[0], fig_file
    )

    # Barcodes: Get counts of initial baseline calls
    crispr_barcode_gene_df = category_counts(
        cell_spot_df,
        gene_cols,
        barcode_cols,
        spot_score_cols,
        spot_parent_cols,
        guide=True,
    )

    # Genes: Get counts of initial baseline calls
    cell_barcode_gene_df = category_counts(
        cell_spot_df,
        gene_cols,
        barcode_cols,
        spot_score_cols,
        spot_parent_cols,
        guide=False,
    )

    # Assign Cell Quality scores based on gene and barcode assignments
    crispr_barcode_gene_df = cell_quality.assign_cell_quality(
        count_df=crispr_barcode_gene_df,
        parent_cols=spot_parent_cols,
        score_col=spot_score_cols[0],
    ).assign(ImageNumber=image_number, site=site)

    cell_barcode_gene_df = cell_quality.assign_cell_quality(
        count_df=cell_barcode_gene_df,
        parent_cols=spot_parent_cols,
        score_col=spot_score_cols[0],
    ).assign(ImageNumber=image_number, site=site)

    # Table 1 - Full Cell and Gene Category with Scores
    out_file = pathlib.PurePath(
        f"{output_dir}/full_cell_category_scores_by_guide.tsv.gz"
    )
    crispr_barcode_gene_df.to_csv(out_file, sep="\t", index=False, compression="gzip")

    # Table 2 - Full Cell and CRISPR Guide Quality Category with Scores
    num_unique_guides = len(
        crispr_barcode_gene_df.loc[:, barcode_cols].squeeze().unique()
    )
    out_file = pathlib.PurePath(f"{output_dir}/full_cell_category_scores.tsv.gz")
    cell_barcode_gene_df.to_csv(out_file, sep="\t", index=False, compression="gzip")

    # Table 3 - Cell Category Summary
    cell_quality_summary_df = cell_quality.summarize_cell_quality_counts(
        quality_df=crispr_barcode_gene_df, parent_cols=spot_parent_cols
    ).assign(ImageNumber=image_number, site=site)

    out_file = pathlib.PurePath(f"{output_dir}/cell_category_summary_count.tsv")
    cell_quality_summary_df.to_csv(out_file, sep="\t", index=False)

    # Table 4 - Gene by cell category counts
    num_unique_genes = len(cell_barcode_gene_df.loc[:, gene_cols].squeeze().unique())

    gene_category_count_df = cell_quality.summarize_perturbation_quality_counts(
        quality_df=crispr_barcode_gene_df,
        parent_cols=spot_parent_cols,
        group_cols=gene_cols,
    ).assign(ImageNumber=image_number, site=site)

    out_file = pathlib.PurePath(f"{output_dir}/gene_by_cell_category_summary_count.tsv")
    gene_category_count_df.to_csv(out_file, sep="\t", index=False)

    # Table 5 - Guide by cell category counts
    guide_category_count_df = cell_quality.summarize_perturbation_quality_counts(
        quality_df=crispr_barcode_gene_df,
        parent_cols=spot_parent_cols,
        group_cols=gene_cols + barcode_cols,
        guide=True,
    ).assign(ImageNumber=image_number, site=site)

    out_file = pathlib.PurePath(
        f"{output_dir}/guide_by_cell_category_summary_count.tsv"
    )
    gene_category_count_df.to_csv(out_file, sep="\t", index=False)

    passed_gene_df = (
        gene_category_count_df.query("Cell_Class in @cell_filter")
        .groupby(gene_cols)["Cell_Count_Per_Gene"]
        .sum()
        .reset_index()
        .sort_values(by="Cell_Count_Per_Gene", ascending=False)
        .reset_index(drop=True)
    )

    passed_gene_df.loc[:, gene_cols] = pd.Categorical(
        passed_gene_df.loc[:, gene_cols].squeeze(),
        categories=passed_gene_df.loc[:, gene_cols].squeeze()
    )

    # Number of non-targetting controls
    num_nt = passed_gene_df.query(
        f"{gene_cols[0]} in @control_barcodes"
    ).Cell_Count_Per_Gene.values[0]

    # Table 6: Complete Site Summary
    descriptive_results = {
        "image_number": image_number,
        "num_unassigned_spots": num_unassigned_spots,
        "num_assigned_spots": num_assigned_spots,
        "num_unique_genes": num_unique_genes,
        "num_unique_guides": num_unique_guides,
        "num_assigned_cells": num_assigned_cells,
        "number_nontarget_controls_good_cells": num_nt,
    }

    output_file = pathlib.PurePath(f"{output_dir}/site_stats.tsv")
    descriptive_results = pd.DataFrame(descriptive_results, index=[site])
    descriptive_results.to_csv(output_file, sep="\t", index=True)
