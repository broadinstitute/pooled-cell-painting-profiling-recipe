import os
import sys
import pathlib
import argparse
import warnings
import logging
import pandas as pd

from pycytominer import feature_select

sys.path.append("config")
from utils import parse_command_args, process_configuration, get_split_aware_site_info

recipe_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(recipe_path, "scripts"))
from io_utils import read_csvs_with_chunksize

logfolder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "logs")
if not os.path.isdir(logfolder):
    os.mkdir(logfolder)
logging.basicConfig(
    filename=os.path.join(logfolder, "3.feature-select.log"),
    level=logging.INFO,
)

args = parse_command_args()

batch_id = args.batch_id
options_config_file = args.options_config_file
experiment_config_file = args.experiment_config_file
split_step = args.split_step

config = process_configuration(
    batch_id,
    step="profile--feature_select",
    options_config=options_config_file,
    experiment_config=experiment_config_file,
)

# Extract config arguments
split_info = config["experiment"]["split"][split_step]
perform = config["options"]["profile"]["feature_select"]["perform"]

# check if this step should be performed
if not perform:
    sys.exit("Config file set to perform=False, not performing {}".format(__file__))

ignore_files = config["options"]["core"]["ignore_files"]
float_format = config["options"]["core"]["float_format"]
compression = config["options"]["core"]["compression"]

input_spotdir = config["directories"]["preprocess"]["spots"]
single_cell_input_dir = config["directories"]["profile"]["single_cell"]
single_cell_file = config["files"]["single_file_only_output_file"]
feature_select_input_dir = config["directories"]["profile"]["profiles"]
feature_select_input_files = config["files"]["normalize_files"]
feature_select_output_files = config["files"]["feature_select_files"]

sc_config = config["options"]["profile"]["single_cell"]
singlecell_from_single_file = sc_config["output_one_single_cell_file_only"]

feature_select_args = config["options"]["profile"]["feature_select"]
feature_select_operations = feature_select_args["operations"]
feature_select_levels = feature_select_args["levels"]
feature_select_drop_samples = feature_select_args["use_samples"]
feature_select_features = feature_select_args["features"]
feature_select_nacutoff = feature_select_args["na_cutoff"]
feature_select_corr_threshold = feature_select_args["corr_threshold"]
force = feature_select_args["force_overwrite"]

print("Starting 3.feature-select.")
logging.info("Starting 3.feature-select.")

sites = [x.name for x in input_spotdir.iterdir() if x.name not in ignore_files]
site_info_dict = get_split_aware_site_info(
    config["experiment"], sites, split_info, separator="___"
)

for data_split_site in site_info_dict:
    for data_level in feature_select_levels:
        if data_level == "single_cell":
            if not singlecell_from_single_file:
                warnings.warn(
                    "Feature select operation is not enabled for site-specific single cell files. Skipping."
                )
                logging.warning(
                    "Feature select operation is not enabled for site-specific single cell files. Skipping."
                )
                continue

        file_to_feature_select = feature_select_input_files[data_level]
        file_to_feature_select = pathlib.Path(
            file_to_feature_select.parents[0],
            file_to_feature_select.name.replace(
                ".csv.gz", f"_{data_split_site}.csv.gz"
            ),
        )

        print(
            f"Now performing feature selection for {data_level}...with operations: {feature_select_operations} for split {data_split_site}"
        )
        logging.info(
            f"Performing feature selection for {data_level} with operations: {feature_select_operations} for split {data_split_site}"
        )

        output_file = feature_select_output_files[data_level]
        output_file = pathlib.Path(
            feature_select_output_files[data_level].parents[0],
            output_file.name.replace(".csv.gz", f"_{data_split_site}.csv.gz"),
        )
        df = read_csvs_with_chunksize(file_to_feature_select)

        feature_select(
            profiles=df,
            features=feature_select_features,
            samples=feature_select_drop_samples,
            operation=feature_select_operations,
            na_cutoff=feature_select_nacutoff,
            corr_threshold=feature_select_corr_threshold,
            output_file=output_file,
            compression_options=compression,
            float_format=float_format,
        )
print("Finished 3.feature-select.")
logging.info("Finished 3.feature-select.")
