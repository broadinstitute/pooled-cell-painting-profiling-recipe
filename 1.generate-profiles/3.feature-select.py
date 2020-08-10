import os
import sys
import pathlib
import argparse
import warnings
import pandas as pd

from pycytominer import feature_select

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
float_format = core_args["float_format"]
compression = core_args["compression"]

normalize_args = config["normalize"]
feature_select_args = config["feature_select"]

singlecell_from_single_file = core_args["output_one_single_cell_file_only"]
feature_select_operations = feature_select_args["operations"]
feature_select_levels = feature_select_args["levels"]
feature_select_drop_samples = feature_select_args["drop_samples"]
feature_select_features = feature_select_args["features"]
feature_select_nacutoff = feature_select_args["na_cutoff"]
feature_select_corr_threshold = feature_select_args["corr_threshold"]
feature_select_input_files = normalize_args["normalize_output_files"]
feature_select_output_files = feature_select_args["feature_select_output_files"]

for data_level in feature_select_levels:
    if data_level == "single_cell":
        if not singlecell_from_single_file:
            warnings.warn(
                "Feature select operation is not enabled for site-specific single cell files. Skipping."
            )
            continue

    input_file = feature_select_input_files[data_level]
    output_file = feature_select_output_files[data_level]

    print(
        f"Now performing feature selection for {data_level}...with operations: {feature_select_operations}"
    )

    feature_select(
        profiles=input_file,
        features=feature_select_features,
        samples=feature_select_drop_samples,
        operation=feature_select_operations,
        na_cutoff=feature_select_nacutoff,
        corr_threshold=feature_select_corr_threshold,
        output_file=output_file,
        compression=compression,
        float_format=float_format,
    )
