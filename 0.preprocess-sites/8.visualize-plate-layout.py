import os
import sys
import pathlib
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
#import numpy as np

sys.path.append(os.path.join("..", "scripts"))
from config_utils import process_config_file
from cell_quality_utils import CellQuality

parser = argparse.ArgumentParser()
parser.add_argument(
    "--config_file",
    help="configuration yaml file for preprocessing pipeline",
    default="site_processing_config.yaml",
)
args = parser.parse_args()
config_file = args.config_file

config = process_config_file(config_file)

# Defines the sections of the config file
core_args = config["core"]
summ_cell_args = config["summarize-cells"]

# Defines the variables set in the config file
batch = core_args["batch"]
quality_func = core_args["categorize_cell_quality"]
ignore_files = core_args["ignore_files"]
acquisition_side = core_args["acquisition_side"]

summ_cells_results_basedir = summ_cell_args["output_resultsdir"]
summ_cells_figures_basedir = summ_cell_args["output_figuresdir"]
output_figuresdir = pathlib.Path(summ_cells_figures_basedir, batch)
cell_category_order = summ_cell_args["cell_category_order"]

cell_count_file = pathlib.Path(summ_cells_results_basedir, batch, "cells", "cell_count.tsv")
cell_count_df = pd.read_csv(cell_count_file, sep="\t")

figures_output = pathlib.Path(summ_cells_figures_basedir, batch)
os.makedirs(figures_output, exist_ok=True)

# creates x, y coordinates for plotting per-plate views.
# assumes image numbering starts in upper left corner and proceeds down
x_seq = 1
y_seq = 1
build_seq = []
final_order = []
for i in range(1, ((acquisition_side * acquisition_side) +1)):
    if len(final_order) == (acquisition_side * acquisition_side):
        break
    else:
        build_seq.append([x_seq, y_seq])
        if len(build_seq) == acquisition_side:
            x_seq = x_seq +1
        if y_seq == acquisition_side:
            y_seq = 0
            build_seq = build_seq[::-1]
            final_order.extend(build_seq)
            build_seq = []
        y_seq = y_seq + 1
sites_list = [*range(1, (acquisition_side * acquisition_side) + 1)]

# uses sites_list in case there are fewer analyzed sites than acquired sites
loc_df = pd.DataFrame(final_order).rename(columns={0:'x_loc', 1:'y_loc'})
loc_df["Site"] = sites_list

# create total_cell_count
cell_count_bysite_df = cell_count_df.groupby("site")["cell_count"].sum().reset_index().rename(columns={'cell_count':'total_cell_count'})
# add total_cell_count to cell_count_df
cell_count_df = cell_count_df.merge(cell_count_bysite_df, on="site")
#add in x, y coordinates for plotting
cell_count_df = cell_count_df.merge(loc_df, left_on="Site", right_on="Site")

# Plot total number of cells per well
cell_count_totalcells_df = cell_count_df.groupby(["x_loc", "y_loc", "Well", "Site"])["total_cell_count"].mean().reset_index()
os.makedirs(output_figuresdir, exist_ok=True)

for well in cell_count_totalcells_df["Well"].unique():
    f, ax = plt.subplots(figsize=(12, 6))
    by_well_df = cell_count_totalcells_df.loc[cell_count_totalcells_df["Well"]==well]
    by_well_df_pivot = by_well_df.pivot(index='y_loc', columns='x_loc', values='total_cell_count')
    labels = by_well_df.pivot(index='y_loc', columns='x_loc', values='Site')
    plt.text(0,0, f"Well {well}", fontsize = 12, color='Black')
    heatmap = sns.heatmap(by_well_df_pivot, annot=labels, fmt = '', ax=ax, xticklabels=False, yticklabels=False)
    fig = heatmap.get_figure()
    output_file = pathlib.Path(output_figuresdir, f"plate_layout_cels_count_per_well_{well}.png")
    fig.savefig(output_file)

# Plot cell number by category per well
f, ax = plt.subplots(len(cell_category_order), len(cell_count_df["Well"].unique()), figsize=(12, 12))
for well in cell_count_df["Well"].unique():
    cell_count_df_slice = cell_count_df.loc[cell_count_df["Well"]==well]
    for quality in cell_category_order:
        all_quality_df = cell_count_df_slice.loc[cell_count_df_slice["Cell_Quality"]==quality]
        single_quality_df = all_quality_df.pivot(index='y_loc', columns='x_loc', values='cell_count')
        labels = all_quality_df.pivot(index='y_loc', columns='x_loc', values='Site')
        sns.heatmap(single_quality_df, ax=ax, fmt = '', xticklabels=False, yticklabels=False)
