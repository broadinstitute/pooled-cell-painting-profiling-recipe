import os
import sys
import pathlib
import argparse
import warnings
import logging
import traceback
import pandas as pd

from pycytominer import aggregate
from pycytominer.cyto_utils import output

sys.path.append("config")
from utils import parse_command_args, process_configuration, get_split_aware_site_info

recipe_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(recipe_path, "scripts"))
from io_utils import read_csvs_with_chunksize

# Configure logging
logfolder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "logs")
if not os.path.isdir(logfolder):
    os.mkdir(logfolder)
logging.basicConfig(
    filename=os.path.join(logfolder, "1.aggregate.log"), level=logging.INFO,
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

config, incomplete_sites, errored_sites = process_configuration(
    batch_id,
    step="profile--aggregate",
    options_config=options_config_file,
    experiment_config=experiment_config_file,
)
logging.info(f"Config used:{config}")
logging.info(f"Skipped incomplete sites during config processing: {incomplete_sites}")
logging.info(f"Skipped errored sites during config processing: {errored_sites}")

# Extract config arguments
split_info = config["experiment"]["split"][split_step]
perform = config["options"]["profile"]["aggregate"]["perform"]

# check if this step should be performed
if not perform:
    sys.exit("Config file set to perform=False, not performing {}".format(__file__))

ignore_files = config["options"]["core"]["ignore_files"]
float_format = config["options"]["core"]["float_format"]
compression = config["options"]["core"]["compression"]

input_spotdir = config["directories"]["preprocess"]["spots"]
single_cell_output_dir = config["directories"]["profile"]["single_cell"]
aggregate_output_dir = config["directories"]["profile"]["profiles"]

single_cell_file = config["files"]["single_file_only_output_file"]
single_cell_site_files = config["files"]["single_cell_site_files"]
aggregate_output_files = config["files"]["aggregate_files"]

sc_config = config["options"]["profile"]["single_cell"]
aggregate_from_single_file = sc_config["output_one_single_cell_file_only"]

aggregate_args = config["options"]["profile"]["aggregate"]
aggregate_operation = aggregate_args["operation"]
aggregate_features = aggregate_args["features"]
aggregate_levels = aggregate_args["levels"]

force = aggregate_args["force_overwrite"]

print("Starting 1.aggregate.")
logging.info(f"Started 1.aggregate.")

sites = [x.name for x in input_spotdir.iterdir() if x.name not in ignore_files]
site_info_dict = get_split_aware_site_info(
    config["experiment"], sites, split_info, separator="___"
)

for data_split_site in site_info_dict:
    # Define a dataset specific file
    single_cell_dataset_file = pathlib.Path(
        single_cell_output_dir,
        single_cell_file.name.replace(".csv.gz", f"_{data_split_site}.csv.gz"),
    )
    # Input argument flow control
    if aggregate_from_single_file:
        assert (
            single_cell_dataset_file.exists()
        ), "Error! The single cell file does not exist! Check 0.merge-single-cells.py"

    # Load single cell data
    if aggregate_from_single_file:
        print(f"Loading one single cell file: {single_cell_dataset_file}")
        single_cell_df = read_csvs_with_chunksize(single_cell_dataset_file, sep=",")
        logging.info(f"Loaded one single cell file: {single_cell_dataset_file}")
    else:
        sites = site_info_dict[data_split_site]
        print(f"Now loading data from {len(sites)} sites")
        logging.info(f"Loading data from {len(sites)} sites")
        single_cell_df = []
        for site in sites:
            site_file = single_cell_site_files[site]
            if site_file.exists():
                site_df = read_csvs_with_chunksize(site_file, sep=",")
                single_cell_df.append(site_df)
            else:
                warnings.warn(
                    f"{site_file} does not exist. There must have been an error in processing"
                )
                logging.warning(f"{site_file} does not exist.")

        single_cell_df = pd.concat(single_cell_df, axis="rows").reset_index(drop=True)

    # Perform the aggregation based on the defined levels and columns
    aggregate_output_dir.mkdir(parents=True, exist_ok=True)
    for aggregate_level, aggregate_columns in aggregate_levels.items():
        aggregate_output_file = aggregate_output_files[aggregate_level]

        print(
            f"Now aggregating by {aggregate_level}...with operation: {aggregate_operation}"
        )
        logging.info(
            f"Aggregating by {aggregate_level}...with operation: {aggregate_operation}"
        )

        aggregate_df = aggregate(
            population_df=single_cell_df,
            strata=aggregate_columns,
            features=aggregate_features,
            operation=aggregate_operation,
        )

        # Define a dataset specific file
        aggregate_dataset_file = pathlib.Path(
            aggregate_output_dir,
            aggregate_output_file.name.replace(".csv.gz", f"_{data_split_site}.csv.gz"),
        )

        output(
            aggregate_df,
            output_filename=aggregate_dataset_file,
            compression_options=compression,
            float_format=float_format,
        )
print("Finished 1.aggregate.")
logging.info(f"Finished 1.aggregate.")
