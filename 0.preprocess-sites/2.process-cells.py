
"""
2.process-cells.py

To view and change parameters used in this file, see `process-cells` in the config.

The following script will process cell painting metadata information for each site.

The metadata include:

1. Cell Counts per Quality Measure
2. Metadata Information per Cell

All sites are processed independently and results are saved in site-specific folders
(in <OUTPUT_BASEDIR>/<BATCH>/paint set in site_processing_config)

1.process-spots must be run before running this script.
"""

import os
import sys
import pathlib
import argparse
import pandas as pd

sys.path.append(os.path.join("..", "scripts"))
from config_utils import process_config_file
from cell_quality_utils import CellQuality
from paint_utils import (
    load_single_cell_compartment_csv,
    merge_single_cell_compartments
)

parser = argparse.ArgumentParser()
parser.add_argument(
    "--config_file",
    help="configuration yaml file for preprocessing pipeline",
    default="site_processing_config.yaml",
)
args = parser.parse_args()
config_file = args.config_file

config = process_config_file(config_file)

#Defines the sections of the config file
core_args = config["core"]
cell_args = config["process-cells"]
prefilter_args = config["prefilter"]
spot_args = config["process-spots"]

#Defines the variables set in the config file
batch = core_args["batch"]
batch_dir = core_args["batch_dir"]
project_dir = core_args["project_dir"]
quality_func = core_args["categorize_cell_quality"]
ignore_files = core_args["ignore_files"]
id_cols = core_args["id_cols"]
compartments = core_args["compartments"]
parent_cols = core_args["parent_cols"]

prefilter_dir = prefilter_args["output_basedir"]
prefilter_file = prefilter_args["prefilter_file"]

foci_dir = spot_args["output_spotdir"]

cell_sort_col = cell_args["sort_col"]
output_basedir = cell_args["output_basedir"]
merge_info = cell_args["merge_columns"]
metadata_merge_left = cell_args["metadata_merge_columns"]["left"]
metadata_merge_right = cell_args["metadata_merge_columns"]["right"]
metadata_merge_cells_left = cell_args["metadata_merge_columns"]["cells_left"]

cell_quality = CellQuality(quality_func)
cell_category_dict = cell_quality.define_cell_quality()
cell_category_df = pd.DataFrame(cell_category_dict, index=["Cell_Class"])

# Enables feature filtering by loading the Cell Painting feature file.
# 0.prefilter-features.py must be run first
try:
    all_feature_df = pd.read_csv(prefilter_file, sep="\t")
except FileNotFoundError:
    raise FileNotFoundError(
        "Error",
        f"{prefilter_file} not found.  ",
        "Perform 0.prefilter-features.py prefilter before continuing..."
    )

sites = [x for x in os.listdir(batch_dir) if x not in ignore_files]

for site in sites:
    try:
        print(f"Now processing {site}...")
        output_folder = pathlib.Path(output_basedir, batch, "paint", site)
        os.makedirs(output_folder, exist_ok=True)

        meta_output_file = pathlib.Path(output_folder, f"metadata_{site}.tsv.gz")
        count_output_file = pathlib.Path(output_folder,f"cell_counts_{site}.tsv")

        compartment_dir = pathlib.Path(batch_dir, site)

        #Make the compartment_csvs dictionary used to merge dfs
        compartment_csvs = {}
        for compartment in compartments:
            try:
                metadata_cols = parent_cols[compartment.lower()] + id_cols
            except KeyError:
                metadata_cols = id_cols
            compartment_csvs[compartment] = load_single_cell_compartment_csv(
                compartment_dir, compartment, metadata_cols
            )

        #Load .tsv created in 1.process-spots
        foci_file = pathlib.Path(foci_dir, site, "full_cell_category_scores_by_guide.tsv.gz")
        try:
            foci_df = pd.read_csv(foci_file, sep="\t")
        except FileNotFoundError:
            print("""Run 1.process-spots before continuing. \n
                  If it has been run, check that output_basedir is set correctly in
                  the process-spots section of the config""")

        #Relabel columns in foci_df to start with Metadata_Foci_
        foci_df.columns = ["Metadata_Foci_{}".format(x) for x in foci_df.columns]

    except FileNotFoundError:
        print(f"{site} data not found")
        continue


    # Merge compartment csvs
    sc_merged_df = merge_single_cell_compartments(compartment_csvs, merge_info, id_cols)
    sc_merged_df = sc_merged_df.sort_values(by=sort_col)

    # Add metadata
    sc_merged_with_metadata_df = foci_df.merge(
        sc_merged_df,
        left_on=metadata_merge_left,
        right_on=metadata_merge_right,
        how="right",
    )

    metadata_df = sc_merged_with_metadata_df.loc[
        :, sc_merged_with_metadata_df.columns.str.contains("Metadata")
    ].drop_duplicates()

    metadata_df.loc[
    metadata_df.Metadata_Foci_Cell_Category.isna(), "Metadata_Foci_Cell_Category"
    ] = 5
    metadata_df.Metadata_Foci_Cell_Category = metadata_df.Metadata_Foci_Cell_Category.astype(
        int
    )
    metadata_df.loc[metadata_df.Metadata_Foci_site.isna(), "Metadata_Foci_site"] = site

    metadata_df = (
    metadata_df.merge(
        cell_category_df,
        left_on=metadata_merge_cells_left,
        right_index=True,
        how="left",
    )
    .sort_values(by=sort_col)
    .reset_index(drop=True)
    .reindex(metadata_col_order, axis="columns")
    )

    #Save files
    metadata_df.to_csv(meta_output_file, sep="\t", index=False)
    cell_count_df.to_csv(count_output_file, sep="\t", index=False)
