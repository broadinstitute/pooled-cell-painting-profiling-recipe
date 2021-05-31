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

sys.path.append("config")
from utils import parse_command_args, process_configuration

recipe_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(recipe_path, "scripts"))
from paint_utils import load_single_cell_compartment_csv, merge_single_cell_compartments
from cell_quality_utils import CellQuality
from profile_utils import sanitize_gene_col

args = parse_command_args()

plate_id = args.plate_id
options_config_file = args.options_config_file
experiment_config_file = args.experiment_config_file

config = process_configuration(
    plate_id,
    step="profile--single_cell",
    options_config=options_config_file,
    experiment_config=experiment_config_file,
)

# Extract config arguments
control_barcodes = config["experiment"]["control_barcode_ids"]

compartments = config["options"]["core"]["compartments"]
cell_filter = config["options"]["core"]["cell_quality"]["cell_filter"]
float_format = config["options"]["core"]["float_format"]
compression = config["options"]["core"]["compression"]
ignore_files = config["options"]["core"]["ignore_files"]
id_cols = config["options"]["core"]["cell_id_cols"]
parent_col_info = config["options"]["core"]["cell_match_cols"]

input_platedir = config["directories"]["input_data_dir"]
single_cell_output_dir = config["directories"]["profile"]["single_cell"]
input_spotdir = config["directories"]["preprocess"]["spots"]
input_paintdir = config["directories"]["preprocess"]["paint"]

prefilter_file = config["files"]["prefilter_file"]
single_file_only_output_file = config["files"]["single_file_only_output_file"]

sc_config = config["options"]["profile"]["single_cell"]
prefilter_features = sc_config["prefilter_features"]
sanitize_genes = sc_config["sanitize_gene_col"]
cell_quality_col = sc_config["cell_quality_column"]
merge_info = sc_config["merge_columns"]
single_file_only = sc_config["output_one_single_cell_file_only"]
force = sc_config["force_overwrite"]
perform = sc_config["perform"]

gene_col = config["options"]["profile"]["aggregate"]["levels"]["gene"]

# check if this step should be performed
if not perform:
    sys.exit("Config file set to perform=False, not performing {}".format(__file__))

# Forced overwrite can be achieved in one of two ways.
# The command line overrides the config file, check here if it is provided
if not force:
    force = args.force

# Check if single cell file already exists, and warn user about no effect
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
sites = [x.name for x in input_spotdir.iterdir() if x.name not in ignore_files]

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
    site_metadata_dir = pathlib.Path(input_paintdir, site)
    site_compartment_dir = pathlib.Path(input_platedir, site)

    # Load cell metadata after cell quality determined in 0.preprocess-sites
    metadata_file = pathlib.Path(site_metadata_dir, f"metadata_{site}.tsv.gz")
    metadata_df = pd.read_csv(metadata_file, sep="\t").query(
        f"{cell_quality_col} in @cell_filter"
    )

    if sanitize_genes:
        metadata_df = sanitize_gene_col(metadata_df, gene_col, control_barcodes)

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
            f"Not all compartments are present in site: {site}\nCheck CellProfiler output path: {site_compartment_dir}. Skipping this site."
        )
        continue

    # Merge single cell compartments together
    sc_merged_df = merge_single_cell_compartments(compartment_csvs, merge_info, id_cols)
    sc_merged_df = sc_merged_df.assign(Metadata_Foci_site=site).reindex(
        ["Metadata_Foci_site"] + all_feature_df.feature_name.tolist(), axis="columns"
    )

    # Merge single cell profiles and metadata info containing cell assignments
    # A left merge ensures that we retain only cells that pass the quality threshold
    common_cols = list(set(metadata_df.columns).intersection(set(sc_merged_df.columns)))
    sc_merged_df = metadata_df.merge(
        sc_merged_df, on=common_cols, how="left"
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
