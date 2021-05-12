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
import pandas as pd

sys.path.append("config")
from utils import parse_command_args, process_configuration

recipe_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(recipe_path, "scripts"))
from cell_quality_utils import CellQuality
from paint_utils import load_single_cell_compartment_csv, merge_single_cell_compartments
from io_utils import check_if_write

args = parse_command_args()

plate_id = args.plate_id
options_config_file = args.options_config_file
experiment_config_file = args.experiment_config_file

config = process_configuration(
    plate_id,
    options_config=options_config_file,
    experiment_config=experiment_config_file,
)

# Define variables set in the config file
ignore_files = config["options"]["core"]["ignore_files"]
id_cols = config["options"]["core"]["cell_id_cols"]
compartments = config["options"]["core"]["compartments"]
parent_cols = config["options"]["core"]["cell_match_cols"]
quality_func = config["options"]["core"]["cell_quality"]["categorize_cell_quality"]
quality_col = config["options"]["core"]["cell_quality"]["cell_quality_column"]
quality_idx = config["options"]["core"]["cell_quality"]["cell_quality_index"]

prefilter_file = config["files"]["prefilter_file"]
input_image_file = config["files"]["image_file"]

input_platedir = config["directories"]["input_data_dir"]
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
    all_feature_df = pd.read_csv(prefilter_file, sep="\t").query("not prefilter_column")
except FileNotFoundError:
    raise FileNotFoundError(
        "Error",
        f"{prefilter_file} not found.  ",
        "Perform 0.prefilter-features.py prefilter before continuing...",
    )

# Load image metadata summary file to extract out important metadata indicators
# 1.process-spots.py must be run first
try:
    image_df = pd.read_csv(input_image_file, sep="\t")
except FileNotFoundError:
    raise FileNotFoundError(
        "Error",
        f"{input_image_file} not found. ",
        "Perform 1.process.spots.py before continuing...",
    )

sites = [x for x in os.listdir(foci_dir) if x not in ignore_files]

for site in sites:
    # Extract out image metadata information for the specific site
    image_subset_df = image_df.query("Metadata_site == @site")

    well = image_subset_df.loc[:, image_cols["well"]].squeeze()
    site_location = image_subset_df.loc[:, image_cols["site"]].squeeze()
    plate = image_subset_df.loc[:, image_cols["plate"]].squeeze()

    try:
        print(f"Now processing cells for {site}...")
        compartment_dir = pathlib.Path(input_platedir, site)

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
            foci_df = pd.read_csv(foci_file, sep="\t")
        except FileNotFoundError:
            print(
                """Run 1.process-spots before continuing. \n
                  If it has been run, check that output_paintdir is set correctly in
                  the process-spots section of the config"""
            )
            continue

        # Relabel columns in foci_df to start with Metadata_Foci_
        foci_df.columns = [
            f"Metadata_Foci_{x}" if not x.startswith("Metadata_Foci") else x
            for x in foci_df.columns
        ]

    except FileNotFoundError:
        print(f"Compartment data for {site} not found")
        continue

    # Merge compartment csvs. Each row is a separate cell
    sc_merged_df = merge_single_cell_compartments(compartment_csvs, merge_info, id_cols)
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
            cell_category_df,
            left_on=quality_idx,
            right_index=True,
            how="left",
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
        )
    )

    output_folder = pathlib.Path(output_paintdir, site)
    if output_folder.exists():
        if force:
            warnings.warn("Output files likely exist, now overwriting...")
        else:
            warnings.warn("Output files likely exist. If they do, NOT overwriting...")
    os.makedirs(output_folder, exist_ok=True)

    meta_output_file = pathlib.Path(output_folder, f"metadata_{site}.tsv.gz")
    count_output_file = pathlib.Path(output_folder, f"cell_counts_{site}.tsv")

    # Save files
    if check_if_write(meta_output_file, force):
        metadata_df.to_csv(meta_output_file, sep="\t", index=False)
    if check_if_write(count_output_file, force):
        cell_count_df.to_csv(count_output_file, sep="\t", index_label="Cell_Quality")

print("All sites complete.")
