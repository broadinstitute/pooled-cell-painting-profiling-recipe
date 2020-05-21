import os
import sys
import pathlib
import argparse
import warnings
import pandas as pd
import yaml

import plotnine as gg
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append("../scripts")
from config_utils import process_config_file
from paint_utils import load_single_cell_compartment_csv, merge_single_cell_compartments
from cell_quality_utils import CellQuality

parser = argparse.ArgumentParser()
parser.add_argument(
    "--config_file",
    help="configuration yaml file for the profiling pipeline",
    default="profiling_config.yaml",
)
parser.add_argument(
    "--single_file_only",
    help="add flag to output all single cell profiles merged into one output file",
    action="store_true",
)
parser.add_argument(
    "--force", help="force overwriting of single cell data", action="store_true"
)
args = parser.parse_args()
config_file = args.config_file
single_file_only = args.single_file_only
force = args.force

config = process_config_file(config_file)

# Extract config arguments
master_args = config["master_config"]
core_args = config["core"]
single_cell_args = config["single_cell"]

project = master_args["project_tag"]
batch = core_args["batch"]
batch_dir = core_args["batch_dir"]
compartments = core_args["compartments"]
parent_col_info = core_args["parent_cols"]
id_cols = core_args["id_cols"]
ignore_files = core_args["ignore_files"]
float_format = core_args["float_format"]
compression = core_args["compression"]

id_cols = core_args["id_cols"]
spot_parent_cols = core_args["parent_cols"]["spots"]

prefilter_features = single_cell_args["prefilter_features"]
prefilter_file = single_cell_args["prefilter_file"]
filter_cell_quality = single_cell_args["filter_cell_quality"]
cell_quality_col = single_cell_args["cell_quality_column"]
spot_batch_dir = single_cell_args["spot_metadata_dir"]
paint_metadata_dir = single_cell_args["paint_metadata_dir"]
merge_info = single_cell_args["merge_columns"]
single_cell_output_dir = single_cell_args["single_cell_output_dir"]
single_file_only_output_file = single_cell_args["single_file_only_output_file"]

if single_file_only:
    if single_file_only_output_file.exists():
        if not force:
            warnings.warn(
                "Combined single cell file exists. Use '--force' to overwrite."
            )

# Load preselected features
all_feature_df = pd.read_csv(prefilter_file, sep="\t")

if prefilter_features:
    all_feature_df = all_feature_df.query("not prefilter_column")

# Pull out all sites that were measured
sites = [x.name for x in spot_batch_dir.iterdir() if x.name not in ignore_files]

sc_df = []
for site in sites:
    # Define single cell output directory and files
    site_output_dir = pathlib.Path(single_cell_output_dir, site)
    site_output_dir.mkdir(parents=True, exist_ok=True)
    sc_output_file = pathlib.Path(site_output_dir, f"{site}_single_cell.csv.gz")

    # Define options based on input flags
    if single_file_only:
        print(f"Building single file; combining single cells from site: {site}...")
    else:
        # If the output file already exists, only overwrite if --force is provided
        if sc_output_file.exists():
            if not force:
                print(
                    f"Skipping reprocessing single cells for site: {site}... use --force to overwrite"
                )
                continue
            else:
                print(f"Now overwriting single cells for site: {site}...")
        else:
            print(f"Now processing single cells for site: {site}...")

    # Point to appropriate directories
    site_metadata_dir = pathlib.Path(paint_metadata_dir, site)
    site_compartment_dir = pathlib.Path(batch_dir, site)

    # Load cell metadata after cell quality determined in 0.preprocess-sites
    metadata_file = pathlib.Path(site_metadata_dir, f"metadata_{site}.tsv.gz")
    metadata_df = pd.read_csv(metadata_file, sep="\t").query(
        f"{cell_quality_col} in @filter_cell_quality"
    )

    # Load csv files for prespecified compartments
    compartment_csvs = {}
    for compartment in compartments:
        try:
            metadata_cols = parent_col_info[compartment.lower()] + id_cols
        except KeyError:
            metadata_cols = id_cols
        try:
            compartment_csvs[compartment] = load_single_cell_compartment_csv(
                site_compartment_dir, compartment, metadata_cols
            )
        except FileNotFoundError:
            continue

    if len(compartment_csvs) != len(compartments):
        warnings.warn(
            f"Not all compartments are present in site: {site}\nCheck CellProfiler output path: {site_compartment_dir}"
        )
        continue

    # Merge single cell compartments together
    sc_merged_df = merge_single_cell_compartments(compartment_csvs, merge_info, id_cols)
    sc_merged_df = sc_merged_df.assign(Metadata_Foci_site=site).reindex(
        ["Metadata_Foci_site"] + all_feature_df.feature_name.tolist(), axis="columns"
    )

    # Merge single cell profiles and metadata
    sc_merged_df = metadata_df.merge(
        sc_merged_df, on=merge_info["metadata_linking_columns"], how="left"
    ).reset_index(drop=True)

    if single_file_only:
        sc_df.append(sc_merged_df)
    else:
        sc_merged_df.to_csv(
            sc_output_file,
            sep=",",
            index=False,
            compression=compression,
            float_format=float_format,
        )

if single_file_only:
    if single_file_only_output_file.exists():
        assert force, "Not outputting combined single cell file, --force not provided!"

    (pd.concat(sc_df, axis="rows").reset_index(drop=True)).to_csv(
        single_file_only_output_file,
        sep=",",
        index=False,
        compression=compression,
        float_format=float_format,
    )
