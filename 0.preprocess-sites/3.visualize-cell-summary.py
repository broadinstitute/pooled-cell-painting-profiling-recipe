import os
import sys
import pathlib
import argparse
import logging
import traceback
import warnings
import pandas as pd
import plotnine as gg
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append("config")
from utils import parse_command_args, process_configuration, get_split_aware_site_info

recipe_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(recipe_path, "scripts"))
from cell_quality_utils import CellQuality
from io_utils import check_if_write, read_csvs_with_chunksize

# Configure logging
logfolder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "logs")
if not os.path.isdir(logfolder):
    os.mkdir(logfolder)
logging.basicConfig(
    filename=os.path.join(logfolder, "3.visualize-cell-summary.log"),
    level=logging.INFO,
)


def handle_excepthook(exc_type, exc_value, exc_traceback):
    logging.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))
    traceback_details = "\n".join(traceback.extract_tb(exc_traceback).format())
    print(f"Uncaught Exception: {traceback_details}")


sys.excepthook = handle_excepthook

# Configure experiment
args = parse_command_args()
logging.info(f"Args used:{args}")

batch_id = args.batch_id
options_config_file = args.options_config_file
experiment_config_file = args.experiment_config_file
split_step = args.split_step

config = process_configuration(
    batch_id,
    step="preprocess--summarize-cells",
    options_config=options_config_file,
    experiment_config=experiment_config_file,
)
logging.info(f"Config used:{config}")

# Define variables set in the config file
control_barcodes = config["experiment"]["control_barcode_ids"]
split_info = config["experiment"]["split"][split_step]

ignore_files = config["options"]["core"]["ignore_files"]
cell_filter = config["options"]["core"]["cell_quality"]["cell_filter"]
quality_func = config["options"]["core"]["cell_quality"]["categorize_cell_quality"]
quality_col = config["options"]["core"]["cell_quality"]["cell_quality_column"]
quality_idx = config["options"]["core"]["cell_quality"]["cell_quality_index"]
cell_category_order = config["options"]["core"]["cell_quality"]["cell_category_order"]
cell_category_colors = config["options"]["core"]["cell_quality"]["cell_category_colors"]

spot_config = config["options"]["preprocess"]["process-spots"]
barcode_cols = spot_config["barcode_cols"]
gene_cols = spot_config["gene_cols"]
spot_score_cols = spot_config["spot_score_cols"]
spot_score_count_cols = ["Metadata_Foci_" + col + "_count" for col in spot_score_cols]
spot_score_mean_cols = ["Metadata_Foci_" + col + "_mean" for col in spot_score_cols]

input_spotdir = config["directories"]["preprocess"]["spots"]
input_paintdir = config["directories"]["preprocess"]["paint"]
output_resultsdir = config["directories"]["preprocess"]["results"]
output_figuresdir = config["directories"]["preprocess"]["figures"]

cell_count_output_file = config["files"]["cell_count_file"]
total_cell_count_file = config["files"]["total_cell_count_file"]

force = config["options"]["preprocess"]["summarize-cells"]["force_overwrite"]
perform = config["options"]["preprocess"]["summarize-cells"]["perform"]

# Perform the pipeline
cell_quality = CellQuality(
    quality_func, category_class_name=quality_col, category_col_index=quality_idx
)
cell_category_dict = cell_quality.define_cell_quality()
empty_cell_category = len(cell_category_dict) + 1
cell_category_dict[empty_cell_category] = "Empty"
cell_category_df = pd.DataFrame(cell_category_dict, index=[quality_col])
cell_category_list = list(cell_category_dict.values())

# Check if this step should be performed
if not perform:
    sys.exit("Config file set to perform=False, not performing {}".format(__file__))

# Forced overwrite can be achieved in one of two ways.
# The command line overrides the config file, check here if it is provided
if not force:
    force = args.force

print("Starting 3.visualize-cell-summary.")
logging.info(f"Starting 3.visualize-cell-summary.")
# Pull out site info and split into distinct datasets based on experiment config
sites = [x.name for x in input_paintdir.iterdir() if x.name not in ignore_files]
print(f"Summarizing {len(sites)} sites in batch: {batch_id}.")
logging.info(f"Summarizing {len(sites)} sites in batch: {batch_id}.")

site_info_dict = get_split_aware_site_info(
    config["experiment"], sites, split_info, separator="___"
)

# Read and Merge Data
cell_quality_list = []
site_stat_list = []
pert_counts_list = []
for data_split_site in site_info_dict:
    split_sites = site_info_dict[data_split_site]

    for site in split_sites:
        # Aggregates cell quality by site into single list
        cell_count_file = pathlib.Path(
            f"{input_paintdir}/{site}/cell_counts_{site}.tsv"
        )
        cell_quality_list.append(read_csvs_with_chunksize(cell_count_file, sep="\t"))

        # Aggregates site summary stats into a single list
        site_stat_file = pathlib.Path(input_spotdir, site, f"site_stats.tsv")
        site_stat_list.append(read_csvs_with_chunksize(site_stat_file, sep="\t"))

        # Aggregates perturbation counts by site into a single list
        pert_count_file = pathlib.Path(
            input_spotdir, site, f"cell_perturbation_category_summary_counts.tsv"
        )
        pert_counts_list.append(read_csvs_with_chunksize(pert_count_file, sep="\t"))

# Creates dataframe from cell quality list
cell_count_df = pd.concat(cell_quality_list, axis="rows").reset_index(drop=True)

# Assigns the Cell_Quality column to the category datatype and sets categories
cell_count_df.loc[:, "Cell_Quality"] = pd.Categorical(
    cell_count_df.Cell_Quality, categories=cell_category_order
)

# Assigns the site column to the category datatype
cell_count_df.loc[:, "site"] = pd.Categorical(
    cell_count_df.site,
    categories=(
        cell_count_df.groupby("site")["cell_count"]
        .sum()
        .sort_values(ascending=False)
        .index.tolist()
    ),
)

if check_if_write(cell_count_output_file, force, throw_warning=True):
    cell_count_df.to_csv(cell_count_output_file, sep="\t", index=False)

# Graph: Cell count with all wells in same graph
cell_count_gg = (
    gg.ggplot(cell_count_df, gg.aes(x="site", y="cell_count"))
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
    cell_count_gg.save(output_file, dpi=300, width=12, height=7, verbose=False)

# Same graph as above, separated by well.
cell_count_gg_parsed = (
    gg.ggplot(cell_count_df, gg.aes(x="site", y="cell_count"))
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
    cell_count_gg_parsed.save(output_file, dpi=300, width=12, height=7, verbose=False)

# Same graph as the two above, separated by dataset splits.
cell_count_gg_dataset_parsed = (
    gg.ggplot(cell_count_df, gg.aes(x="site", y="cell_count"))
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
    + gg.facet_wrap("~Metadata_dataset_split", drop=False, scales="free_x")
)

output_file = pathlib.Path(
    output_figuresdir, "all_cellpainting_cellquality_across_sites_by_datasetsplit.png"
)
if check_if_write(output_file, force, throw_warning=True):
    cell_count_gg_dataset_parsed.save(
        output_file, dpi=300, width=12, height=7, verbose=False
    )

#  Total cells in each quality category
all_count_df = pd.DataFrame(
    cell_count_df.groupby(["Cell_Quality", "Metadata_dataset_split"])[
        "cell_count"
    ].sum()
).reset_index()
if check_if_write(total_cell_count_file, force, throw_warning=True):
    all_count_df.to_csv(total_cell_count_file, sep="\t", index=False)

# Graph: Total cells in each quality category
all_cells = all_count_df.cell_count.sum()

total_cell_count_gg = (
    gg.ggplot(all_count_df, gg.aes(x="Cell_Quality", y="cell_count"))
    + gg.geom_bar(gg.aes(fill="Cell_Quality"), stat="identity")
    + gg.geom_text(gg.aes(label="cell_count"))
    + gg.theme_bw()
    + gg.theme(
        axis_text_x=gg.element_text(rotation=90, size=9),
        strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
    )
    + gg.xlab("")
    + gg.ylab("Cell Count")
    + gg.ggtitle(f"{all_cells} Total Cells")
    + gg.facet_wrap("~Metadata_dataset_split", drop=False, scales="free_x")
    + gg.scale_fill_manual(
        name="Cell Quality", labels=cell_category_order, values=cell_category_colors,
    )
)

output_file = pathlib.Path(output_figuresdir, "total_cell_count.png")
if check_if_write(output_file, force, throw_warning=True):
    total_cell_count_gg.save(output_file, dpi=300, width=7, height=8, verbose=False)

# Total cell number by well
all_well_count_df = pd.DataFrame(
    cell_count_df.groupby(["Cell_Quality", "well"])["cell_count"].sum()
).reset_index()

# Graph: Total cell number by well
total_cell_well_count_gg = (
    gg.ggplot(all_well_count_df, gg.aes(x="well", y="cell_count"))
    + gg.geom_bar(gg.aes(fill="Cell_Quality"), stat="identity")
    + gg.geom_text(gg.aes(label="cell_count"))
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
    id_vars=["site", "Metadata_dataset_split", "image_number"],
    value_vars=["num_unique_genes", "num_unique_guides"],
    value_name="pert_count",
    var_name="pert_class",
)

# Assigns the site column to the category datatype
num_unique_pert_df.loc[:, "site"] = pd.Categorical(
    num_unique_pert_df.site,
    categories=(
        site_stat_df.groupby("site")["num_assigned_cells"]
        .sum()
        .sort_values(ascending=False)
        .index.tolist()
    ),
)

unique_pert_count_gg = (
    gg.ggplot(num_unique_pert_df, gg.aes(x="site", y="pert_count"))
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
pert_count_df = pd.concat(pert_counts_list, axis="rows").reset_index()

# Output a full count of perturbations per site
output_file = pathlib.Path(
    output_resultsdir, "complete_perturbation_count_per_site.tsv.gz"
)

if check_if_write(output_file, force, throw_warning=True):
    pert_count_df.to_csv(output_file, index=False, sep="\t")

# Summarize counts further in preparation for plotting
pert_count_df = (
    pert_count_df.reset_index(drop=True)
    .groupby(gene_cols + barcode_cols + [quality_col])["Cell_Count_Per_Guide"]
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

# Assigns a category datatype
pert_count_df.loc[:, gene_cols] = pd.Categorical(
    pert_count_df.loc[:, gene_cols].squeeze(),
    categories=(
        pert_count_df.groupby(gene_cols)["Cell_Count_Per_Guide"]
        .sum()
        .sort_values(ascending=False)
        .index.tolist()
    ),
)

output_file = pathlib.Path(
    output_resultsdir, "all_cellpainting_unique_perturbations_coverage.tsv"
)
if check_if_write(output_file, force, throw_warning=True):
    pert_count_df.to_csv(output_file, index=False, sep="\t")

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

# Drop controls from the visualization
pert_count_no_control_df = pert_count_df.loc[
    ~pert_count_df.loc[:, gene_cols].isin(control_barcodes).squeeze(),
].reset_index(drop=True)

guide_count_no_control_gg = (
    gg.ggplot(
        pert_count_no_control_df, gg.aes(x=gene_cols[0], y="Cell_Count_Per_Guide")
    )
    + gg.geom_bar(gg.aes(fill="factor(barcode_id)"), stat="identity")
    + gg.theme_bw()
    + gg.theme(axis_text_y=gg.element_text(size=5))
    + gg.xlab("Genes")
    + gg.ylab("Coverage (Cell Count)")
    + gg.coord_flip()
    + gg.scale_fill_discrete(name="Guide")
)

output_file = pathlib.Path(output_figuresdir, "perturbation_coverage_no_controls.png")
if check_if_write(output_file, force, throw_warning=True):
    guide_count_no_control_gg.save(
        output_file, dpi=300, width=7, height=20, verbose=False
    )

print(f"There are a total of {all_cells} cells in {batch_id}")
logging.info(f"There are a total of {all_cells} cells in {batch_id}")
print("Finished 3.visualize-cell-summary.")
logging.info(f"Finished 3.visualize-cell-summary.")
