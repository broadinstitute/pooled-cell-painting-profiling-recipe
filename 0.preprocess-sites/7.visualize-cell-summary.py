import os
import sys
import pathlib
import argparse
import pandas as pd
import plotnine as gg
import matplotlib.pyplot as plt
import seaborn as sns

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
spot_args = config["process-spots"]
cell_args = config["process-cells"]
summ_cell_args = config["summarize-cells"]

# Defines the variables set in the config file
batch = core_args["batch"]
quality_func = core_args["categorize_cell_quality"]

barcode_cols = spot_args["barcode_cols"]
barcode_cols = ["Metadata_Foci_" + col for col in barcode_cols]
gene_cols = spot_args["gene_cols"]
gene_cols = ["Metadata_Foci_" + col for col in gene_cols]
spot_score_cols = spot_args["spot_score_cols"]
spot_score_count_cols = ["Metadata_Foci_" + col + "_count" for col in spot_score_cols]
spot_score_mean_cols = ["Metadata_Foci_" + col + "_mean" for col in spot_score_cols]

input_basedir = cell_args["output_basedir"]
metadata_foci_col = cell_args["metadata_merge_columns"]["cell_quality_col"]
cell_cols = cell_args["metadata_merge_columns"]["cell_cols"]
cell_quality_col = cell_args["metadata_merge_columns"]["cell_quality_col"]
foci_site_col = cell_args["foci_site_col"]

output_resultsdir = summ_cell_args["output_resultsdir"]
output_resultsdir = pathlib.Path(output_resultsdir, batch)
output_figuresdir = summ_cell_args["output_figuresdir"]
output_figuresdir = pathlib.Path(output_figuresdir, batch)
cell_category_order = summ_cell_args["cell_category_order"]
cell_category_colors = summ_cell_args["cell_category_colors"]

cell_quality = CellQuality(quality_func)
cell_category_dict = cell_quality.define_cell_quality()
empty_cell_category = len(cell_category_dict) + 1
cell_category_dict[empty_cell_category] = "Empty"
cell_category_df = pd.DataFrame(cell_category_dict, index=["Cell_Class"])
cell_category_list = list(cell_category_dict.values())

# Read and Merge Data
cell_quality_list = []
metadata_list = []
metadata_col_list = (
    ["Cell_Class"]
    + cell_cols
    + barcode_cols
    + gene_cols
    + spot_score_count_cols
    + spot_score_mean_cols
    + [cell_quality_col, foci_site_col]
)

input_dir = pathlib.Path(input_basedir, batch, "paint")
sites = [x for x in os.listdir(input_dir) if x not in ignore_files]
print(f"There are {len(sites)} sites.")

for site in sites:

    cell_count_file = pathlib.Path(input_dir, site, f"cell_counts_{site}.tsv")
    metadata_file = pathlib.Path(input_dir, site, f"metadata_{site}.tsv.gz")

    # Aggregates cell quality by site into single list
    cell_quality_list.append(pd.read_csv(cell_count_file, sep="\t"))

    # Aggregates metadata by site into a single list
    metadata_df = (
        pd.read_csv(metadata_file, sep="\t")
        .loc[:, metadata_col_list]
        .reset_index(drop=True)
    )

    metadata_list.append(metadata_df)

# Creates dataframe from cell quality list
cell_count_df = (
    pd.concat(cell_quality_list, axis="rows")
    .rename(columns={"site": "site_full"})
    .reset_index(drop=True)
)

# Assigns the Cell_Quality column to the category datatype and sets categories
cell_count_df.loc[:, "Cell_Quality"] = pd.Categorical(
    cell_count_df.Cell_Quality, categories=cell_category_order
)

# Assigns the site_full column to the category datatype
cell_count_df.loc[:, "site_full"] = pd.Categorical(
    cell_count_df.site_full,
    categories=(
        cell_count_df.groupby("site_full")["cell_count"]
        .sum()
        .sort_values(ascending=False)
        .index.tolist()
    ),
)

cell_count_df = cell_count_df.assign(
    Plate=[x[0] for x in cell_count_df.site_full.str.split("-")],
    Well=[x[1] for x in cell_count_df.site_full.str.split("-")],
    Site=[x[2] for x in cell_count_df.site_full.str.split("-")],
)

output_folder = pathlib.Path(output_resultsdir, "cells")
os.makedirs(output_folder, exist_ok=True)
output_file = pathlib.Path(output_folder, "cell_count.tsv")
cell_count_df.to_csv(output_file, sep="\t", index=False)


# Graph: Cell count with all wells in same graph
cell_count_gg = (
    gg.ggplot(cell_count_df, gg.aes(x="site_full", y="cell_count"))
    + gg.geom_bar(gg.aes(fill="Cell_Quality"), stat="identity")
    + gg.theme_bw()
    + gg.theme(axis_text_x=gg.element_text(rotation=90, size=5))
    + gg.xlab("Sites")
    + gg.ylab("Cell Count")
    + gg.scale_fill_manual(
        name="Cell Quality", labels=cell_category_list, values=cell_category_colors
    )
)

os.makedirs(output_figuresdir, exist_ok=True)
output_file = pathlib.Path(
    output_figuresdir, "all_cellpainting_cellquality_across_sites.png"
)
cell_count_gg.save(output_file, dpi=300, width=10, height=7, verbose=False)

# Same graph as above, separated by Well.
cell_count_gg_parsed = (
    gg.ggplot(cell_count_df, gg.aes(x="site_full", y="cell_count"))
    + gg.geom_bar(gg.aes(fill="Cell_Quality"), stat="identity")
    + gg.theme_bw()
    + gg.theme(
        axis_text_x=gg.element_text(rotation=90, size=5),
        strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
    )
    + gg.xlab("Sites")
    + gg.ylab("Cell Count")
    + gg.scale_fill_manual(
        name="Cell Quality", labels=cell_category_order, values=cell_category_colors
    )
    + gg.facet_wrap("~Well", drop=False, scales="free_x")
)

output_file = pathlib.Path(
    output_figuresdir, "all_cellpainting_cellquality_across_sites_by_well.png"
)
cell_count_gg_parsed.save(output_file, dpi=300, width=10, height=7, verbose=False)

#  Total cells in each quality category
all_count_df = pd.DataFrame(
    cell_count_df.groupby("Cell_Quality")["cell_count"].sum()
).reset_index()
output_folder = pathlib.Path(output_resultsdir, "cells")
output_file = pathlib.Path(output_folder, "total_cell_count.tsv")
all_count_df.to_csv(output_file, sep="\t", index=False)

# Graph: Total cells in each quality category
all_cells = all_count_df.cell_count.sum()

total_cell_count_gg = (
    gg.ggplot(all_count_df, gg.aes(x="Cell_Quality", y="cell_count"))
    + gg.geom_bar(gg.aes(fill="Cell_Quality"), stat="identity")
    + gg.theme_bw()
    + gg.theme(axis_text_x=gg.element_text(rotation=90, size=9))
    + gg.xlab("")
    + gg.ylab("Cell Count")
    + gg.ggtitle(f"{all_cells} Total Cells")
    + gg.scale_fill_manual(
        name="Cell Quality", labels=cell_category_order, values=cell_category_colors,
    )
)

output_file = pathlib.Path(output_figuresdir, "total_cell_count.png")
total_cell_count_gg.save(output_file, dpi=300, width=5, height=6, verbose=False)

total_cell_count_gg

print(f"There are a total of {all_cells} cells in {batch}")

# Total cell number by well
all_well_count_df = pd.DataFrame(
    cell_count_df.groupby(["Cell_Quality", "Well"])["cell_count"].sum()
).reset_index()

# Graph: Total cell number by well
total_cell_well_count_gg = (
    gg.ggplot(all_well_count_df, gg.aes(x="Well", y="cell_count"))
    + gg.geom_bar(gg.aes(fill="Cell_Quality"), stat="identity")
    + gg.theme_bw()
    + gg.theme(axis_text_x=gg.element_text(rotation=90, size=9))
    + gg.xlab("")
    + gg.ylab("Cell Count")
    + gg.facet_wrap("~Cell_Quality")
    + gg.scale_fill_manual(
        name="Cell Quality", labels=cell_category_order, values=cell_category_colors,
    )
    + gg.theme(strip_background=gg.element_rect(colour="black", fill="#fdfff4"))
)

output_file = pathlib.Path(output_figuresdir, "total_cell_count_by_well.png")
total_cell_well_count_gg.save(output_file, dpi=400, width=6, height=5, verbose=False)
