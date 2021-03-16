import os
import sys
import pathlib
import argparse
import warnings
import pandas as pd

from pycytominer import aggregate
from pycytominer.cyto_utils import output

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
perform = config["options"]["profile"]["aggregate"]["perform"]

# check if this step should be performed
if not perform:
    sys.exit("Config file set to perform=False, not performing {}".format(__file__))

ignore_files = config["options"]["core"]["ignore_files"]
float_format = config["options"]["core"]["float_format"]
compression = config["options"]["core"]["compression"]

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

# Input argument flow control
if aggregate_from_single_file:
    assert (
        single_cell_file.exists()
    ), "Error! The single cell file does not exist! Check 0.merge-single-cells.py"

# Load single cell data
if aggregate_from_single_file:
    print(f"Loading one single cell file: {single_cell_file}")
    single_cell_df = pd.read_csv(single_cell_file, sep=",")
else:
    sites = list(single_cell_site_files)
    print(f"Now loading data from {len(sites)} sites")
    single_cell_df = []
    for site in sites:
        site_file = single_cell_site_files[site]
        if site_file.exists():
            site_df = pd.read_csv(site_file, sep=",")
            single_cell_df.append(site_df)
        else:
            warnings.warn(
                f"{site_file} does not exist. There must have been an error in processing"
            )

    single_cell_df = pd.concat(single_cell_df, axis="rows").reset_index(drop=True)

# Perform the aggregation based on the defined levels and columns
aggregate_output_dir.mkdir(parents=True, exist_ok=True)
for aggregate_level, aggregate_columns in aggregate_levels.items():
    aggregate_output_file = aggregate_output_files[aggregate_level]

    print(
        f"Now aggregating by {aggregate_level}...with operation: {aggregate_operation}"
    )

    aggregate_df = aggregate(
        population_df=single_cell_df,
        strata=aggregate_columns,
        features=aggregate_features,
        operation=aggregate_operation,
    )

    output(
        aggregate_df,
        output_filename=aggregate_output_file,
        compression=compression,
        float_format=float_format,
    )
