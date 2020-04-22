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
trash_files = core_args["trash_files"]

quality_func = spot_args["cell_quality"]
output_basedir = spot_args["output_basedir"]

id_cols = core_args["id_cols"]
spot_parent_cols = core_args["parent_cols"]["spots"]
barcode_cols = spot_args["barcode_cols"]
gene_cols = spot_args["gene_cols"]
location_cols = spot_args["location_cols"]
spot_score_cols = spot_args["spot_score_cols"]
foci_cols = spot_args["foci_cols"]

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
        print("{} data not found".format(site))

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
        print("{} data not aligned between foci files")

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
    crispr_barcode_gene_count_df = category_counts(
        cell_spot_df,
        gene_cols,
        barcode_cols,
        spot_score_cols,
        spot_parent_cols,
        guide=True,
    )

    # Genes: Get counts of initial baseline calls
    cell_barcode_gene_count_df = category_counts(
        cell_spot_df,
        gene_cols,
        barcode_cols,
        spot_score_cols,
        spot_parent_cols,
        guide=False,
    )
