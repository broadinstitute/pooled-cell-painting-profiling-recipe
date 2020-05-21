import os
import sys
import pathlib
import argparse
import warnings
import pandas as pd

from pycytominer import aggregate
from pycytominer.cyto_utils import output

sys.path.append("../scripts")
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
ignore_files = core_args["ignore_files"]
float_format = core_args["float_format"]
compression = core_args["compression"]

single_cell_args = config["single_cell"]
single_cell_output_dir = single_cell_args["single_cell_output_dir"]
single_cell_file = single_cell_args["single_file_only_output_file"]

aggregate_args = config["aggregate"]
aggregate_from_site = aggregate_args["from_site_files"]
aggregate_from_single_file = aggregate_args["from_single_file"]
aggregate_output_dir = aggregate_args["aggregate_output_dir"]
aggregate_operation = aggregate_args["operation"]
aggregate_features = aggregate_args["features"]
aggregate_levels = aggregate_args["levels"]

# Input argument flow control
if aggregate_from_single_file:
    assert single_cell_file.exists(), "Error! The master single cell file does not exist! Did you mean to set `from_single_file` to False?"
else:
    assert aggregate_from_site, "You must set one aggregate source to true in the config (either `from_site_files` or `from_single_file`)"

if aggregate_from_single_file and aggregate_from_site:
    warnings.warn(
        "Both sources set to true. Defaulting to aggregate `from_single_file`"
    )

# Load single cell data
if aggregate_from_single_file:
    single_cell_df = pd.read_csv(single_cell_file, sep=",")
else:
    sites = [x for x in single_cell_output_dir.iterdir() if x.name not in ignore_files]
    single_cell_df = []
    for site in sites:
        if site.is_file():
            continue
        site_file = pathlib.Path(site, f"{site.name}_single_cell.csv.gz")
        if site_file.exists():
            site_df = pd.read_csv(site_file, sep=",")
            single_cell_df.append(site_df)
        else:
            warnings.warn(f"{site_file} does not exist. There must have been an error in processing")

    single_cell_df = pd.concat(single_cell_df, axis="rows").reset_index(drop=True)

# Perform the aggregation based on the defined levels and columns
aggregate_output_dir.mkdir(parents=True, exist_ok=True)
for aggregate_level, aggregate_columns in aggregate_levels.items():
    aggregate_output_file = pathlib.Path(aggregate_output_dir, f"{batch}_{aggregate_level}.csv.gz")

    print(f"Now aggregating by {aggregate_level}...with operation: {aggregate_operation}")

    aggregate_df = aggregate(
        population_df=single_cell_df,
        strata=aggregate_columns,
        features=aggregate_features,
        operation=aggregate_operation
    )

    output(
        aggregate_df,
        output_filename=aggregate_output_file,
        compression=compression,
        float_format=float_format
    )
