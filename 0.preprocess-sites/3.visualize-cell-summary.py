import os
import sys
import pathlib
import argparse
import warnings
import pandas as pd
import plotnine as gg
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append("config")
from config_utils import process_config_file

recipe_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(recipe_path, "scripts"))
from cell_quality_utils import CellQuality
from arg_utils import parse_command_args
from io_utils import check_if_write

args = parse_command_args(config_file="site_processing_config.yaml")
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
ignore_files = core_args["ignore_files"]
control_barcodes = core_args["control_barcodes"]

barcode_cols = spot_args["barcode_cols"]
gene_cols = spot_args["gene_cols"]
spot_score_cols = spot_args["spot_score_cols"]
spot_score_count_cols = ["Metadata_Foci_" + col + "_count" for col in spot_score_cols]
spot_score_mean_cols = ["Metadata_Foci_" + col + "_mean" for col in spot_score_cols]
input_paintdir = cell_args["output_paintdir"]
input_spotdir = cell_args["output_spotdir"]

output_resultsdir = summ_cell_args["output_summary_resultsdir"]
output_figuresdir = summ_cell_args["output_summary_figuresdir"]
cell_count_output_file = summ_cell_args["cell_count_file"]
total_cell_count_file = summ_cell_args["total_cell_count_file"]
cell_category_order = summ_cell_args["cell_category_order"]
cell_category_colors = summ_cell_args["cell_category_colors"]

cell_quality = CellQuality(quality_func)
cell_category_dict = cell_quality.define_cell_quality()
empty_cell_category = len(cell_category_dict) + 1
cell_category_dict[empty_cell_category] = "Empty"
cell_category_df = pd.DataFrame(cell_category_dict, index=["Cell_Class"])
cell_category_list = list(cell_category_dict.values())
force = summ_cell_args["force_overwrite"]

# Forced overwrite can be achieved in one of two ways.
# The command line overrides the config file, check here if it is provided
if not force:
    force = args.force

sites = [x.name for x in input_paintdir.iterdir() if x.name not in ignore_files]
print(f"There are {len(sites)} sites.")

# Read and Merge Data
cell_quality_list = []
site_stat_list = []
pert_counts_list = []
for site in sites:
    # Aggregates cell quality by site into single list
    cell_count_file = pathlib.Path(f"{input_paintdir}/{site}/cell_counts_{site}.tsv")
    cell_quality_list.append(pd.read_csv(cell_count_file, sep="\t"))

    # Aggregates site summary stats into a single list
    site_stat_file = pathlib.Path(spot_input_dir, site, f"site_stats.tsv")
    site_stat_list.append(pd.read_csv(site_stat_file, sep="\t"))

    # Aggregates perturbation counts by site into a single list
    pert_count_file = pathlib.Path(
        spot_input_dir, site, f"cell_perturbation_category_summary_counts.tsv"
    )
    pert_counts_list.append(pd.read_csv(pert_count_file, sep="\t"))

# Creates dataframe from cell quality list
cell_count_df = pd.concat(cell_quality_list, axis="rows").reset_index(drop=True)

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

if check_if_write(cell_count_output_file, force, throw_warning=True):
    cell_count_df.to_csv(cell_count_output_file, sep="\t", index=False)

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
if check_if_write(output_file, force, throw_warning=True):
    cell_count_gg.save(output_file, dpi=300, width=10, height=7, verbose=False)

# Same graph as above, separated by well.
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
    + gg.facet_wrap("~well", drop=False, scales="free_x")
)

output_file = pathlib.Path(
    output_figuresdir, "all_cellpainting_cellquality_across_sites_by_well.png"
)
if check_if_write(output_file, force, throw_warning=True):
    cell_count_gg_parsed.save(output_file, dpi=300, width=10, height=7, verbose=False)

#  Total cells in each quality category
all_count_df = pd.DataFrame(
    cell_count_df.groupby("Cell_Quality")["cell_count"].sum()
).reset_index()
if check_if_write(total_cell_count_file, force, throw_warning=True):
    all_count_df.to_csv(total_cell_count_file, sep="\t", index=False)

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
if check_if_write(output_file, force, throw_warning=True):
    total_cell_count_gg.save(output_file, dpi=300, width=5, height=6, verbose=False)

# Total cell number by well
all_well_count_df = pd.DataFrame(
    cell_count_df.groupby(["Cell_Quality", "well"])["cell_count"].sum()
).reset_index()

# Graph: Total cell number by well
total_cell_well_count_gg = (
    gg.ggplot(all_well_count_df, gg.aes(x="well", y="cell_count"))
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
if check_if_write(output_file, force, throw_warning=True):
    total_cell_well_count_gg.save(
        output_file, dpi=400, width=6, height=5, verbose=False
    )

# Process summary statistics per site
site_stat_df = pd.concat(site_stat_list, axis="rows").reset_index(drop=True)

num_unique_pert_df = site_stat_df.melt(
    id_vars=["site_full", "image_number"],
    value_vars=["num_unique_genes", "num_unique_guides"],
    value_name="pert_count",
    var_name="pert_class",
)

# Assigns the site_full column to the category datatype
num_unique_pert_df.loc[:, "site_full"] = pd.Categorical(
    num_unique_pert_df.site_full,
    categories=(
        site_stat_df.groupby("site_full")["num_assigned_cells"]
        .sum()
        .sort_values(ascending=False)
        .index.tolist()
    ),
)

unique_pert_count_gg = (
    gg.ggplot(num_unique_pert_df, gg.aes(x="site_full", y="pert_count"))
    + gg.geom_bar(gg.aes(fill="pert_class"), stat="identity")
    + gg.theme_bw()
    + gg.theme(axis_text_x=gg.element_text(rotation=90, size=5))
    + gg.xlab("Sites")
    + gg.ylab("Perturbation Count")
    + gg.scale_fill_discrete(name="Perturbation Class")
)

output_file = pathlib.Path(
    output_figuresdir, "all_cellpainting_unique_perturbations_across_sites.png"
)
if check_if_write(output_file, force, throw_warning=True):
    unique_pert_count_gg.save(output_file, dpi=300, width=10, height=7, verbose=False)

# Process overall perturbation counts per batch
pert_count_df = (
    pd.concat(pert_counts_list, axis="rows")
    .reset_index()
    .rename(columns={"site": "site_full"})
)
pert_count_df = pert_count_df.loc[
    ~pert_count_df.loc[:, gene_cols].isin(control_barcodes).squeeze(),
].reset_index(drop=True)
pert_count_df = (
    pert_count_df.groupby(gene_cols + barcode_cols)["Cell_Count_Per_Guide"]
    .sum()
    .reset_index()
)

id_group_df = (
    pert_count_df.groupby(gene_cols)[barcode_cols]
    .transform(lambda x: pd.factorize(x)[0])
    .reset_index(drop=True)
)
id_group_df.columns = ["barcode_id"]

pert_count_df = pert_count_df.merge(id_group_df, left_index=True, right_index=True)

# Assigns the site_full column to the category datatype
pert_count_df.loc[:, gene_cols] = pd.Categorical(
    pert_count_df.loc[:, gene_cols].squeeze(),
    categories=(
        pert_count_df.groupby(gene_cols)["Cell_Count_Per_Guide"]
        .sum()
        .sort_values(ascending=False)
        .index.tolist()
    ),
)

guide_count_gg = (
    gg.ggplot(pert_count_df, gg.aes(x=gene_cols[0], y="Cell_Count_Per_Guide"))
    + gg.geom_bar(gg.aes(fill="factor(barcode_id)"), stat="identity")
    + gg.theme_bw()
    + gg.theme(axis_text_y=gg.element_text(size=5))
    + gg.xlab("Genes")
    + gg.ylab("Coverage (Cell Count)")
    + gg.coord_flip()
    + gg.scale_fill_discrete(name="Guide")
)

output_file = pathlib.Path(output_figuresdir, "perturbation_coverage.png")
if check_if_write(output_file, force, throw_warning=True):
    guide_count_gg.save(output_file, dpi=300, width=7, height=20, verbose=False)

print(f"There are a total of {all_cells} cells in {batch}")
