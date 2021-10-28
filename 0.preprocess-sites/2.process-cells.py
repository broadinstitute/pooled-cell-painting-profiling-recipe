"""
2.process-cells.py

To view and change parameters used in this file, see `process-cells` in the config.

The following script will process cell painting metadata information for each site.

The metadata include:

1. Cell Counts per Quality Measure
2. Metadata Information per Cell

All sites are processed independently and results are saved in site-specific folders
(in <OUTPUT_BASEDIR>/<PLATE_ID>/paint set in site_processing_config)

1.process-spots must be run before running this script.
"""

import os
import sys
import pathlib
import warnings
import argparse
import logging
import traceback
import pandas as pd

sys.path.append("config")
from utils import parse_command_args, process_configuration, get_split_aware_site_info

recipe_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(recipe_path, "scripts"))
from cell_quality_utils import CellQuality
from paint_utils import load_single_cell_compartment_csv, merge_single_cell_compartments
from io_utils import check_if_write, read_csvs_with_chunksize

# Configure logging
logfolder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "logs")
if not os.path.isdir(logfolder):
    os.mkdir(logfolder)
logging.basicConfig(
    filename=os.path.join(logfolder, "2.process-cells.log"), level=logging.INFO,
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
    step="preprocess--process-cells",
    options_config=options_config_file,
    experiment_config=experiment_config_file,
)
logging.info(f"Config used:{config}")

# Define variables set in the config file
split_info = config["experiment"]["split"][split_step]

ignore_files = config["options"]["core"]["ignore_files"]
id_cols = config["options"]["core"]["cell_id_cols"]
compartments = config["options"]["core"]["compartments"]
parent_cols = config["options"]["core"]["cell_match_cols"]
quality_func = config["options"]["core"]["cell_quality"]["categorize_cell_quality"]
quality_col = config["options"]["core"]["cell_quality"]["cell_quality_column"]
quality_idx = config["options"]["core"]["cell_quality"]["cell_quality_index"]

prefilter_file = config["files"]["prefilter_file"]
input_image_file = config["files"]["image_file"]

input_batchdir = config["directories"]["input_data_dir"]
foci_dir = config["directories"]["preprocess"]["spots"]
output_paintdir = config["directories"]["preprocess"]["paint"]

image_cols = config["options"]["preprocess"]["process-spots"]["image_cols"]

cell_config = config["options"]["preprocess"]["process-cells"]
cell_sort_col = cell_config["sort_col"]
merge_info = cell_config["merge_columns"]
foci_site_col = cell_config["foci_site_col"]
force = cell_config["force_overwrite"]
perform = cell_config["perform"]
metadata_merge_foci_cols = cell_config["metadata_merge_columns"]["foci_cols"]
metadata_merge_cell_cols = cell_config["metadata_merge_columns"]["cell_cols"]

# check if this step should be performed
if not perform:
    sys.exit("Config file set to perform=False, not performing {}".format(__file__))

# Forced overwrite can be achieved in one of two ways.
# The command line overrides the config file, check here if it is provided
if not force:
    force = args.force

print("Starting 2.process-cells.")
logging.info(f"Starting 2.process-cells.")
cell_quality = CellQuality(
    quality_func, category_class_name=quality_col, category_col_index=quality_idx
)
cell_category_dict = cell_quality.define_cell_quality()
empty_cell_category = len(cell_category_dict) + 1
cell_category_dict[empty_cell_category] = "Empty"
cell_category_df = pd.DataFrame(cell_category_dict, index=[quality_col]).transpose()

# Enables feature filtering by loading the Cell Painting feature file.
# 0.prefilter-features.py must be run first
try:
    all_feature_df = read_csvs_with_chunksize(prefilter_file, sep="\t").query(
        "not prefilter_column"
    )
except FileNotFoundError:
    raise FileNotFoundError(
        "Error",
        f"{prefilter_file} not found.  ",
        "Perform 0.prefilter-features.py prefilter before continuing...",
    )

# Load image metadata summary file to extract out important metadata indicators
# 1.process-spots.py must be run first
try:
    image_df = read_csvs_with_chunksize(input_image_file, sep="\t")
except FileNotFoundError:
    raise FileNotFoundError(
        "Error",
        f"{input_image_file} not found. ",
        "Perform 1.process.spots.py before continuing...",
    )

# Pull out site info and split into distinct datasets based on experiment config
sites = [x for x in os.listdir(foci_dir) if x not in ignore_files]
site_info_dict = get_split_aware_site_info(
    config["experiment"], sites, split_info, separator="___"
)

for data_split_site in site_info_dict:
    split_sites = site_info_dict[data_split_site]

    for site in split_sites:
        # Extract out image metadata information for the specific site
        image_subset_df = image_df.query(f"{image_cols['full_info']} == @site")

        well = image_subset_df.loc[:, image_cols["well"]].squeeze()
        site_location = image_subset_df.loc[:, image_cols["site"]].squeeze()
        plate = image_subset_df.loc[:, image_cols["plate"]].squeeze()

        try:
            print(f"Now processing cells for {site}...")
            logging.info(f"Processing cells for {site}")
            compartment_dir = pathlib.Path(input_batchdir, site)

            # Make the compartment_csvs dictionary used to merge dfs
            compartment_csvs = {}
            for compartment in compartments:
                try:
                    metadata_cols = parent_cols[compartment.lower()] + id_cols
                except KeyError:
                    metadata_cols = id_cols
                compartment_csvs[compartment] = load_single_cell_compartment_csv(
                    compartment_dir, compartment, metadata_cols
                )

            # Load .tsv created in 1.process-spots
            foci_file = pathlib.Path(
                foci_dir, site, "cell_id_barcode_alignment_scores_by_guide.tsv.gz"
            )
            try:
                foci_df = read_csvs_with_chunksize(foci_file, sep="\t")
            except FileNotFoundError:
                print(
                    """Run 1.process-spots before continuing. \n
                      If it has been run, check that output_paintdir is set correctly in
                      the process-spots section of the config"""
                )
                logging.info(
                    f"Didn't find cell_id_barcode_alignment_scores_by_guide.tsv"
                )
                continue

            # Relabel columns in foci_df to start with Metadata_Foci_
            foci_df.columns = [
                f"Metadata_Foci_{x}" if not x.startswith("Metadata_Foci") else x
                for x in foci_df.columns
            ]

        except FileNotFoundError:
            print(f"Compartment data for {site} not found")
            logging.info(f"Didn't find compartment data for {site}")
            continue

        # Merge compartment csvs. Each row is a separate cell
        sc_merged_df = merge_single_cell_compartments(
            compartment_csvs, merge_info, id_cols
        )
        sc_merged_df = sc_merged_df.sort_values(by=cell_sort_col).reindex(
            all_feature_df.feature_name.tolist(), axis="columns"
        )

        # Add metadata to the merged csvs and prefilter to remove sneaky flagged columns
        sc_merged_with_metadata_df = foci_df.merge(
            sc_merged_df,
            left_on=metadata_merge_foci_cols,
            right_on=metadata_merge_cell_cols,
            how="right",
        )

        # Drops all columns that are not Metadata (i.e. profiles)
        metadata_df = sc_merged_with_metadata_df.loc[
            :, sc_merged_with_metadata_df.columns.str.contains("Metadata")
        ].drop_duplicates()

        # Adds a cell quality category to previously uncategorized cells
        metadata_df.loc[:, quality_idx] = metadata_df.loc[:, quality_idx].fillna(
            empty_cell_category
        )

        # Adds the site to the metadata_foci_site column to previously uncategorized cells
        metadata_df[foci_site_col] = metadata_df[foci_site_col].fillna(site)

        # Add the cell quality metadata to the df
        metadata_df = (
            metadata_df.merge(
                cell_category_df, left_on=quality_idx, right_index=True, how="left",
            )
            .sort_values(by=cell_sort_col)
            .drop_duplicates(subset=[cell_sort_col, quality_idx])
            .reset_index(drop=True)
        )

        assert (
            not metadata_df.loc[:, cell_sort_col].duplicated().any()
        ), "Stop! You're counting cells more than once"

        # Create a summary of counts of each cell quality class
        cell_count_df = (
            pd.DataFrame(metadata_df.loc[:, quality_col].value_counts())
            .rename(columns={quality_col: "cell_count"})
            .assign(
                site=site,
                plate=plate,
                well=well,
                site_location=site_location,
                Metadata_dataset_split=data_split_site,
            )
        )

        output_folder = pathlib.Path(output_paintdir, site)
        if output_folder.exists():
            if force:
                warnings.warn("Output files likely exist, now overwriting...")
                logging.warning("Output files likely exist. Overwriting.")
            else:
                warnings.warn(
                    "Output files likely exist. If they do, NOT overwriting..."
                )
                logging.warning("Output files likely exist. NOT Overwriting.")
        os.makedirs(output_folder, exist_ok=True)

        meta_output_file = pathlib.Path(output_folder, f"metadata_{site}.tsv.gz")
        count_output_file = pathlib.Path(output_folder, f"cell_counts_{site}.tsv")

        # Save files
        if check_if_write(meta_output_file, force):
            metadata_df.to_csv(meta_output_file, sep="\t", index=False)
        if check_if_write(count_output_file, force):
            cell_count_df.to_csv(
                count_output_file, sep="\t", index_label="Cell_Quality"
            )

print("Finished 2.process-cells.")
logging.info(f"Finished 2.process-cells.")
