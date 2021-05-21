import os
import sys
import pathlib
import argparse
import warnings
import pandas as pd

from pycytominer import normalize

sys.path.append("config")
from utils import parse_command_args, process_configuration

args = parse_command_args()

plate_id = args.plate_id
options_config_file = args.options_config_file
experiment_config_file = args.experiment_config_file

config = process_configuration(
    plate_id,
    step="profile--normalize",
    options_config=options_config_file,
    experiment_config=experiment_config_file,
)

# Extract config arguments
perform = config["options"]["profile"]["normalize"]["perform"]

# check if this step should be performed
if not perform:
    sys.exit("Config file set to perform=False, not performing {}".format(__file__))

float_format = config["options"]["core"]["float_format"]
compression = config["options"]["core"]["compression"]

normalize_input_files = config["files"]["aggregate_files"]
normalize_output_files = config["files"]["normalize_files"]

sc_config = config["options"]["profile"]["single_cell"]
normalize_singlecell_from_single_file = sc_config["output_one_single_cell_file_only"]

normalize_args = config["options"]["profile"]["normalize"]
normalize_levels = normalize_args["levels"]
normalize_by_samples = normalize_args["by_samples"]
normalize_these_features = normalize_args["features"]
normalize_method = normalize_args["method"]
force = normalize_args["force_overwrite"]

for data_level in normalize_levels:
    if data_level == "single_cell":
        if not normalize_singlecell_from_single_file:
            continue

    file_to_normalize = normalize_input_files[data_level]
    output_file = normalize_output_files[data_level]

    print(f"Now normalizing {data_level}...with operation: {normalize_method}")

    df = pd.read_csv(file_to_normalize)

    normalize(
        profiles=df,
        features=normalize_these_features,
        samples=normalize_by_samples,
        method=normalize_method,
        output_file=output_file,
        compression=compression,
        float_format=float_format,
    )
