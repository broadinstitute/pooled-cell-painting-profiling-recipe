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
import warnings
import argparse
import logging
import traceback
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
    filename=os.path.join(logfolder, "1.process-spots.log"), level=logging.INFO,
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
    step="preprocess--process-spots",
    options_config=options_config_file,
    experiment_config=experiment_config_file,
)
logging.info(f"Config used:{config}")

# Set constants
control_barcodes = config["experiment"]["control_barcode_ids"]
split_info = config["experiment"]["split"][split_step]

id_cols = config["options"]["core"]["cell_id_cols"]
spot_parent_cols = config["options"]["core"]["cell_match_cols"]["spots"]
ignore_files = config["options"]["core"]["ignore_files"]
cell_filter = config["options"]["core"]["cell_quality"]["cell_filter"]
quality_func = config["options"]["core"]["cell_quality"]["categorize_cell_quality"]
quality_col = config["options"]["core"]["cell_quality"]["cell_quality_column"]
quality_idx = config["options"]["core"]["cell_quality"]["cell_quality_index"]

output_image_file = config["files"]["image_file"]
output_spotdir = config["directories"]["preprocess"]["spots"]
input_batchdir = config["directories"]["input_data_dir"]

spot_config = config["options"]["preprocess"]["process-spots"]
image_cols = spot_config["image_cols"]
barcode_cols = spot_config["barcode_cols"]
gene_cols = spot_config["gene_cols"]
location_cols = spot_config["location_cols"]
spot_score_cols = spot_config["spot_score_cols"]
foci_cols = spot_config["foci_cols"]
force = spot_config["force_overwrite"]
perform = spot_config["perform"]
allowed_skips = spot_config["allowed_skips"]

# check if this step should be performed
if not perform:
    sys.exit("Config file set to perform=False, not performing {}".format(__file__))

# Forced overwrite can be achieved in one of two ways.
# The command line overrides the config file, check here if it is provided
if not force:
    force = args.force

print("Starting 1.process-spots.")
logging.info("Started 1.process-spots.")
barcode_foci_cols = id_cols + location_cols + spot_parent_cols
all_foci_cols = list(
    set(
        id_cols + location_cols + foci_cols + spot_score_cols + barcode_cols + gene_cols
    )
)

cell_quality = CellQuality(
    quality_func, category_class_name=quality_col, category_col_index=quality_idx
)
cell_category_dict = cell_quality.define_cell_quality()
cell_category_df = pd.DataFrame(cell_category_dict, index=[quality_col])

# Define site and data set splits
sites = [x.name for x in input_batchdir.iterdir() if x.name not in ignore_files]
num_sites = len(sites)
site_info_dict = get_split_aware_site_info(
    config["experiment"], sites, split_info, separator="___"
)

image_list = []
all_called_barcodes = []
allowed_skip_counter = 0
for data_split_site in site_info_dict:
    split_sites = site_info_dict[data_split_site]
    for site in split_sites:
        if allowed_skips >= allowed_skip_counter:
            print(f"Now processing spots for {site}...part of set {data_split_site}")
            logging.info(f"Now processing spots for {site}...part of set {data_split_site}")

            # Load image metadata per site
            try:
                image_file = pathlib.Path(input_batchdir, site, "Image.csv")
                image_df = read_csvs_with_chunksize(image_file).assign(
                    Metadata_site=site, Metadata_dataset_split=data_split_site
                )
                image_list.append(image_df)

                # Obtain specific metadata info
                well = image_df.loc[:, image_cols["well"]].squeeze()
                plate = image_df.loc[:, image_cols["plate"]].squeeze()
                site_location = image_df.loc[:, image_cols["site"]].squeeze()
            except FileNotFoundError:
                print(f"{site} image metadata does not exist. Skipping...")
                logging.info(f"Skipped {site}. No Image.csv")
                continue
            except:
                print(f"Couldn't parse {site} image metadata. Skipping...")
                logging.warning(f"Couldn't parse {site} image metadata. Skipping...")
                allowed_skip_counter += 1
                print(f"Now at {allowed_skip_counter} sites skipped from errors.")
                logging.warning(f"Now at {allowed_skip_counter} sites skipped from errors.")
                logging.warning(traceback.format_exc())
                continue

            # Load spot data
            try:
                barcode_file = pathlib.Path(input_batchdir, site, "BarcodeFoci.csv")
                barcodefoci_df = read_csvs_with_chunksize(barcode_file)

                foci_file = pathlib.Path(input_batchdir, site, "Foci.csv")
                foci_df = read_csvs_with_chunksize(foci_file)
            except FileNotFoundError:
                print(f"{site} data not found")
                logging.info(f"Skipped {site}. No Foci.csv and/or BarcodeFoci.csv")
                continue
            except:
                print(f"Couldn't parse {site} foci data. Skipping...")
                logging.warning(f"Couldn't parse {site} foci data. Skipping...")
                allowed_skip_counter += 1
                print(f"Now at {allowed_skip_counter} sites skipped from errors.")
                logging.warning(f"Now at {allowed_skip_counter} sites skipped from errors.")
                logging.warning(traceback.format_exc())
                continue

            try:
                image_number = foci_df.ImageNumber.unique()[0]
            except IndexError:
                print(f"{site} does not have any foci")
                logging.info(f"{site} does not have any foci")
                continue

            try:
                # Confirm that image number and object number are aligned
                pd.testing.assert_frame_equal(
                    barcodefoci_df.loc[:, id_cols],
                    foci_df.loc[:, id_cols],
                    check_names=True,
                )

                pd.testing.assert_frame_equal(
                    barcodefoci_df.loc[:, location_cols],
                    foci_df.loc[:, location_cols],
                    check_names=True,
                )
            except AssertionError:
                print(f"{site} data not aligned between foci files")
                logging.info(f"{site} data not aligned between Barcode.csv and Foci.csv")
                continue

            output_dir = pathlib.Path(output_spotdir, site)
            if output_dir.exists():
                if force:
                    warnings.warn("Output files likely exist, now overwriting...")
                    logging.warning(f"{site} output files likely exist. Overwriting.")
                else:
                    warnings.warn(
                        "Output files likely exist. If they do, NOT overwriting..."
                    )
                    logging.warning(f"{site} output files likely exist. NOT Overwriting.")
            output_dir.mkdir(exist_ok=True, parents=True)

            # Merge spot data files
            complete_foci_df = barcodefoci_df.loc[:, barcode_foci_cols].merge(
                foci_df.loc[:, all_foci_cols],
                left_on=id_cols + location_cols,
                right_on=id_cols + location_cols,
                how="inner",
            )

            null_spot_df = complete_foci_df.loc[
                (complete_foci_df.loc[:, spot_parent_cols] == 0).squeeze(), :
            ]
            cell_spot_df = complete_foci_df.loc[
                (complete_foci_df.loc[:, spot_parent_cols] != 0).squeeze(), :
            ]
            num_assigned_cells = len(
                cell_spot_df.loc[:, spot_parent_cols].squeeze().unique()
            )
            num_unassigned_spots = null_spot_df.shape[0]
            num_assigned_spots = cell_spot_df.shape[0]

            # Figure 1 - histogram of barcode counts per cell
            fig_file = pathlib.Path(output_dir, "num_spots_per_cell_histogram.png")
            if check_if_write(fig_file, force):
                spot_counts_per_cell_histogram(cell_spot_df, spot_parent_cols, fig_file)

            # Figure 2 - histogram of barcode scores per spot
            fig_file = pathlib.Path(output_dir, "barcode_scores_per_spot_histogram.png")
            if check_if_write(fig_file, force):
                spot_score_histogram(cell_spot_df, spot_score_cols, fig_file)

            # Figure 3 - Joint plot of relationship of barcode counts per cell and mean score
            fig_file = pathlib.Path(
                output_dir, "per_cell_barcode_count_by_mean_score_jointplot.png"
            )
            if check_if_write(fig_file, force):
                spot_count_score_jointplot(
                    cell_spot_df, spot_parent_cols[0], spot_score_cols[0], fig_file
                )

            # Barcodes: Get counts of initial baseline calls
            crispr_barcode_gene_df = category_counts(
                df=cell_spot_df,
                gene_cols=gene_cols,
                barcode_cols=barcode_cols,
                score_cols=spot_score_cols,
                parent_cols=spot_parent_cols,
                guide=True,
            )

            # Assign Cell Quality scores based on gene and barcode assignments
            crispr_barcode_gene_df = cell_quality.assign_cell_quality(
                count_df=crispr_barcode_gene_df,
                parent_cols=spot_parent_cols,
                score_col=spot_score_cols[0],
            ).assign(
                ImageNumber=image_number,
                site=site,
                plate=plate,
                well=well,
                site_location=site_location,
                Metadata_dataset_split=data_split_site,
            )

            num_unique_guides = len(
                crispr_barcode_gene_df.loc[:, barcode_cols].squeeze().unique()
            )
            num_unique_genes = len(
                crispr_barcode_gene_df.loc[:, gene_cols].squeeze().unique()
            )

            # Table 1 - Full cell and CRISPR guide quality with scores
            out_file = pathlib.Path(
                output_dir, "cell_id_barcode_alignment_scores_by_guide.tsv.gz"
            )
            if check_if_write(out_file, force):
                crispr_barcode_gene_df.to_csv(
                    out_file, sep="\t", index=False, compression="gzip"
                )

            # Table 2 - Cell Category Summary
            cell_quality_summary_df = cell_quality.summarize_cell_quality_counts(
                quality_df=crispr_barcode_gene_df, parent_cols=spot_parent_cols
            ).assign(
                ImageNumber=image_number,
                site=site,
                plate=plate,
                well=well,
                site_location=site_location,
                Metadata_dataset_split=data_split_site,
            )

            out_file = pathlib.Path(output_dir, "cell_category_summary_count.tsv")
            if check_if_write(out_file, force):
                cell_quality_summary_df.to_csv(out_file, sep="\t", index=False)

            # Table 3 - Counting gene and guide by cell category
            gene_category_count_df = cell_quality.summarize_perturbation_quality_counts(
                quality_df=crispr_barcode_gene_df,
                parent_cols=spot_parent_cols,
                group_cols=gene_cols,
            )

            guide_category_count_df = cell_quality.summarize_perturbation_quality_counts(
                quality_df=crispr_barcode_gene_df,
                parent_cols=spot_parent_cols,
                group_cols=gene_cols + barcode_cols,
                guide=True,
            )

            count_merge_cols = list(
                set(gene_category_count_df.columns).intersection(
                    guide_category_count_df.columns
                )
            )

            cell_category_counts_df = (
                guide_category_count_df.merge(
                    gene_category_count_df, on=count_merge_cols, how="left"
                )
                .assign(
                    ImageNumber=image_number,
                    site=site,
                    plate=plate,
                    well=well,
                    site_location=site_location,
                    Metadata_dataset_split=data_split_site,
                )
                .query(f"{quality_col} in @cell_filter")
            )

            out_file = pathlib.Path(
                output_dir, "cell_perturbation_category_summary_counts.tsv"
            )
            if check_if_write(out_file, force):
                cell_category_counts_df.to_csv(out_file, sep="\t", index=False)

            passed_gene_df = (
                gene_category_count_df.groupby(gene_cols)["Cell_Count_Per_Gene"]
                .sum()
                .reset_index()
                .sort_values(by="Cell_Count_Per_Gene", ascending=False)
                .reset_index(drop=True)
            )

            passed_gene_df.loc[:, gene_cols] = pd.Categorical(
                passed_gene_df.loc[:, gene_cols].squeeze(),
                categories=passed_gene_df.loc[:, gene_cols].squeeze(),
            )

            # Number of non-targetting controls
            nt_gene_df = passed_gene_df.query(f"{gene_cols[0]} in @control_barcodes")
            num_nt = nt_gene_df.Cell_Count_Per_Gene.sum()

            # Table 4: Complete Site Summary
            descriptive_results = {
                "image_number": image_number,
                "num_unassigned_spots": num_unassigned_spots,
                "num_assigned_spots": num_assigned_spots,
                "num_unique_genes": num_unique_genes,
                "num_unique_guides": num_unique_guides,
                "num_assigned_cells": num_assigned_cells,
                "number_nontarget_controls_good_cells": num_nt,
            }

            descriptive_results = pd.DataFrame(descriptive_results, index=[0]).assign(
                ImageNumber=image_number,
                site=site,
                plate=plate,
                well=well,
                site_location=site_location,
                Metadata_dataset_split=data_split_site,
            )

            output_file = pathlib.Path(output_dir, "site_stats.tsv")
            if check_if_write(output_file, force):
                descriptive_results.to_csv(output_file, sep="\t", index=False)

            # Append site to barcode count summary
            site_called_barcodes = list(crispr_barcode_gene_df["Barcode_MatchedTo_Barcode"])
            all_called_barcodes += site_called_barcodes
        else:
            print (f"More than {allowed_skips} sites errored and were skipped. Processing stopped.")
            logging.warning(f"More than {allowed_skips} sites errored and were skipped. Processing stopped.")
            sys.exit()

# Create and save the barcode count summary for easy NGS comparison
barcode_count_summary_df = pd.DataFrame()
barcode_count_summary_df["All_Called_Barcodes"] = all_called_barcodes
barcode_count_summary_df = (
    barcode_count_summary_df.groupby(["All_Called_Barcodes"])
    .size()
    .reset_index(name="count")
)
output_file = pathlib.Path(output_spotdir, "Total_Barcode_Calls_Counts.tsv")
barcode_count_summary_df.to_csv(output_file, sep="\t", index=False)

# Save the output image info
image_df = pd.concat(image_list, axis="rows").reset_index(drop=True)
image_df.to_csv(output_image_file, sep="\t", index=False)
print("Finished 1.process-spots.")
logging.info(f"Finished 1.process-spots.")
