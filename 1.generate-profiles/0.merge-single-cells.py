import os
import sys
import pathlib
import argparse
import warnings
import logging
import traceback
import pandas as pd
import yaml

import plotnine as gg
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append("config")
from utils import parse_command_args, process_configuration, get_split_aware_site_info

recipe_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(recipe_path, "scripts"))
from paint_utils import load_single_cell_compartment_csv, merge_single_cell_compartments
from cell_quality_utils import CellQuality
from profile_utils import sanitize_gene_col
from io_utils import read_csvs_with_chunksize

# Configure logging
logfolder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "logs")
if not os.path.isdir(logfolder):
    os.mkdir(logfolder)
logging.basicConfig(
    filename=os.path.join(logfolder, "0.merge-single-cells.log"), level=logging.INFO,
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
    step="profile--single_cell",
    options_config=options_config_file,
    experiment_config=experiment_config_file,
)
logging.info(f"Config used:{config}")

# Extract config arguments
control_barcodes = config["experiment"]["control_barcode_ids"]
split_info = config["experiment"]["split"][split_step]

compartments = config["options"]["core"]["compartments"]
cell_filter = config["options"]["core"]["cell_quality"]["cell_filter"]
float_format = config["options"]["core"]["float_format"]
compression = config["options"]["core"]["compression"]
ignore_files = config["options"]["core"]["ignore_files"]
id_cols = config["options"]["core"]["cell_id_cols"]
parent_col_info = config["options"]["core"]["cell_match_cols"]

input_batchdir = config["directories"]["input_data_dir"]
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
            logging.warning("Combined single cell file already exists.")

print("Starting 0.merge-single-cells.")
logging.info(f"Started 0.merge-single-cells.")

# Load preselected features
all_feature_df = read_csvs_with_chunksize(prefilter_file, sep="\t")

if prefilter_features:
    all_feature_df = all_feature_df.query("not prefilter_column")

# Pull out all sites that were measured
sites = [x.name for x in input_spotdir.iterdir() if x.name not in ignore_files]
site_info_dict = get_split_aware_site_info(
    config["experiment"], sites, split_info, separator="___"
)

for data_split_site in site_info_dict:
    split_sites = site_info_dict[data_split_site]

    sc_df = []
    for site in split_sites:
        # Define single cell output directory and files
        site_output_dir = pathlib.Path(single_cell_output_dir, site)
        site_output_dir.mkdir(parents=True, exist_ok=True)
        sc_output_file = pathlib.Path(site_output_dir, f"{site}_single_cell.csv.gz")

        # Define options based on input flags
        if single_file_only:
            print(
                f"Building single file for dataset {data_split_site}; combining single cells from site: {site}..."
            )
            logging.info(
                f"Building single file for dataset {data_split_site}; combining single cells from site: {site}..."
            )
        else:
            # If the output file already exists, only overwrite if --force is provided
            if sc_output_file.exists():
                if not force:
                    print(
                        f"Skipping reprocessing single cells for site: {site}... use --force to overwrite"
                    )
                    logging.info(f"Skipped reprocessing single cells for site: {site}")
                    continue
                else:
                    print(f"Now overwriting single cells for site: {site}...")
                    logging.info(f"Overwrote single cells for site: {site}")
            else:
                print(f"Now processing single cells for site: {site}...")
                logging.info(f"Processed single cells for site: {site}")

        # Point to appropriate directories
        site_metadata_dir = pathlib.Path(input_paintdir, site)
        site_compartment_dir = pathlib.Path(input_batchdir, site)

        # Load cell metadata after cell quality determined in 0.preprocess-sites
        metadata_file = pathlib.Path(site_metadata_dir, f"metadata_{site}.tsv.gz")
        metadata_df = read_csvs_with_chunksize(metadata_file, sep="\t").query(
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
            logging.warning(
                f"{site} skipped because of missing compartments in {site_compartment_dir}."
            )
            continue

        # Merge single cell compartments together
        sc_merged_df = merge_single_cell_compartments(
            compartment_csvs, merge_info, id_cols
        )
        sc_merged_df = sc_merged_df.assign(Metadata_Foci_site=site).reindex(
            ["Metadata_Foci_site"] + all_feature_df.feature_name.tolist(),
            axis="columns",
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
        # Define a dataset specific file
        single_cell_dataset_file = pathlib.Path(
            single_cell_output_dir,
            single_file_only_output_file.name.replace(
                ".csv.gz", f"_{data_split_site}.csv.gz"
            ),
        )
        if single_cell_dataset_file.exists():
            assert (
                force
            ), "Not outputting combined single cell file, --force not provided!"

        print(f"Now writing single cell file: {single_cell_dataset_file}...")
        logging.info(f"Writing single cell file: {single_cell_dataset_file}...")
        (pd.concat(sc_df, axis="rows").reset_index(drop=True)).to_csv(
            single_cell_dataset_file,
            sep=",",
            index=False,
            compression=compression,
            float_format=float_format,
        )
print("Finished 0.merge-single-cells.")
logging.info(f"Finished 0.merge-single-cells.")
