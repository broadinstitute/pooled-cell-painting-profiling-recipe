import os
import sys
import pathlib
import argparse
import warnings
import pandas as pd

from pycytominer import normalize

sys.path.append("config")
from config_utils import process_config_file

parser = argparse.ArgumentParser()
parser.add_argument(
    "--config_file",
    help="configuration yaml file for the profiling pipeline",
    default="profiling_config.yaml",
)
args = parser.parse_args()
config_file = args.config_file

config = process_config_file(config_file)

# Extract config arguments
core_args = config["core"]
batch = core_args["batch"]
aggregate_args = config["aggregate"]
normalize_args = config["normalize"]

ignore_files = core_args["ignore_files"]
float_format = core_args["float_format"]
compression = core_args["compression"]

normalize_singlecell_from_single_file = core_args["output_one_single_cell_file_only"]
normalize_levels = normalize_args["levels"]
normalize_by_samples = normalize_args["by_samples"]
normalize_these_features = normalize_args["features"]
normalize_method = normalize_args["method"]
normalize_input_files = aggregate_args["aggregate_output_files"]
normalize_output_files = normalize_args["normalize_output_files"]

for data_level in normalize_levels:
    if data_level == "single_cell":
        if not normalize_singlecell_from_single_file:
            continue

    file_to_normalize = normalize_input_files[data_level]
    output_file = normalize_output_files[data_level]

    print(f"Now normalizing {data_level}...with operation: {normalize_method}")

    normalize_df = normalize(
        profiles=file_to_normalize,
        features=normalize_these_features,
        samples=normalize_by_samples,
        method=normalize_method,
        output_file=output_file,
        compression=compression,
        float_format=float_format,
    )
