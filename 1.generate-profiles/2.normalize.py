import os
import sys
import pathlib
import argparse
import warnings
import logging
import traceback
import pandas as pd

from pycytominer import normalize

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
    filename=os.path.join(logfolder, "2.normalize.log"), level=logging.INFO,
)


def handle_excepthook(exc_type, exc_value, exc_traceback):
    logging.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))
    traceback_details = "\n".join(traceback.extract_tb(exc_traceback).format())
    print(f"Uncaught Exception: {traceback_details}")


sys.excepthook = handle_excepthook

# Configure experiment
args = parse_command_args()

batch_id = args.batch_id
options_config_file = args.options_config_file
experiment_config_file = args.experiment_config_file
split_step = args.split_step

config = process_configuration(
    batch_id,
    step="profile--normalize",
    options_config=options_config_file,
    experiment_config=experiment_config_file,
)
logging.info(f"Config used:{config}")

# Extract config arguments
split_info = config["experiment"]["split"][split_step]
perform = config["options"]["profile"]["normalize"]["perform"]

# check if this step should be performed
if not perform:
    sys.exit("Config file set to perform=False, not performing {}".format(__file__))

ignore_files = config["options"]["core"]["ignore_files"]
float_format = config["options"]["core"]["float_format"]
compression = config["options"]["core"]["compression"]

input_spotdir = config["directories"]["preprocess"]["spots"]
normalize_input_dir = config["directories"]["profile"]["profiles"]
single_cell_input_dir = config["directories"]["profile"]["single_cell"]
normalize_input_files = config["files"]["aggregate_files"]
normalize_output_files = config["files"]["normalize_files"]
single_cell_file = config["files"]["single_file_only_output_file"]

sc_config = config["options"]["profile"]["single_cell"]
normalize_singlecell_from_single_file = sc_config["output_one_single_cell_file_only"]

normalize_args = config["options"]["profile"]["normalize"]
normalize_levels = normalize_args["levels"]
normalize_by_samples = normalize_args["by_samples"]
normalize_these_features = normalize_args["features"]
normalize_method = normalize_args["method"]
force = normalize_args["force_overwrite"]

print("Starting 2.normalize.")
logging.info(f"Started 2.normalize.")

sites = [x.name for x in input_spotdir.iterdir() if x.name not in ignore_files]
site_info_dict = get_split_aware_site_info(
    config["experiment"], sites, split_info, separator="___"
)

for data_split_site in site_info_dict:
    for data_level in normalize_levels:
        if data_level == "single_cell":
            if not normalize_singlecell_from_single_file:
                continue
            file_to_normalize = pathlib.Path(
                single_cell_input_dir,
                single_cell_file.name.replace(".csv.gz", f"_{data_split_site}.csv.gz"),
            )
        else:
            file_to_normalize = normalize_input_files[data_level]
            file_to_normalize = pathlib.Path(
                normalize_input_dir,
                file_to_normalize.name.replace(".csv.gz", f"_{data_split_site}.csv.gz"),
            )

        print(
            f"Now normalizing {data_level}...with operation: {normalize_method} for split {data_split_site}"
        )
        logging.info(
            f"Normalizing {data_level}...with operation: {normalize_method} for split {data_split_site}"
        )

        output_file = normalize_output_files[data_level]
        output_file = pathlib.Path(
            normalize_output_files[data_level].parents[0],
            output_file.name.replace(".csv.gz", f"_{data_split_site}.csv.gz"),
        )
        df = read_csvs_with_chunksize(file_to_normalize)

        normalize(
            profiles=df,
            features=normalize_these_features,
            samples=normalize_by_samples,
            method=normalize_method,
            output_file=output_file,
            compression_options=compression,
            float_format=float_format,
        )
print("Finished 2.normalize.")
logging.info("Finished 2.normalize.")
