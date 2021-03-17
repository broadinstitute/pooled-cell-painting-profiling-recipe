import os
import sys
import pathlib
import argparse
import warnings
import pandas as pd

from pycytominer import feature_select

sys.path.append("config")
from utils import parse_command_args, process_configuration

args = parse_command_args()

plate_id = args.plate_id
options_config_file = args.options_config_file
experiment_config_file = args.experiment_config_file

config = process_configuration(
    plate_id,
    options_config=options_config_file,
    experiment_config=experiment_config_file,
)

# Extract config arguments
perform = config["options"]["profile"]["feature_select"]["perform"]
# check if this step should be performed
if not perform:
    sys.exit("Config file set to perform=False, not performing {}".format(__file__))

float_format = config["options"]["core"]["float_format"]
compression = config["options"]["core"]["compression"]

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

    df = pd.read_csv(input_file)

    feature_select(
        profiles=df,
        features=feature_select_features,
        samples=feature_select_drop_samples,
        operation=feature_select_operations,
        na_cutoff=feature_select_nacutoff,
        corr_threshold=feature_select_corr_threshold,
        output_file=output_file,
        compression=compression,
        float_format=float_format,
    )
