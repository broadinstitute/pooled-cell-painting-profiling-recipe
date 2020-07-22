"""
0.prefilter-features.py

Determine which features should be used for building morphology profiles.

Note this is a preselection step.
An additional round of feature selection will occur at a different stage.
"""

import os
import sys
import pathlib
import argparse
import numpy as np
import pandas as pd
from scripts.site_processing_utils import prefilter_features

sys.path.append(os.path.join("..", "scripts"))
from config_utils import process_config_file

parser = argparse.ArgumentParser()
parser.add_argument(
    "--config_file",
    help="configuration yaml file for preprocessing pipeline",
    default="site_processing_config.yaml",
)
parser.add_argument(
    "--force", help="force overwriting of feature data", action="store_true"
)
args = parser.parse_args()
config_file = args.config_file
force = args.force

config = process_config_file(config_file)

# Set constants
main_args = config["main_config"]
core_args = config["core"]
prefilter_args = config["prefilter"]

project = main_args["project_tag"]
batch = core_args["batch"]
compartments = core_args["compartments"]

perform = prefilter_args["perform"]
example_site_dir = prefilter_args["example_site_dir"]
flag_cols = prefilter_args["flag_cols"]
output_file = prefilter_args["prefilter_file"]
force = prefilter_args["force_overwrite"]

# Forced overwrite can be achieved in one of two ways.
# The command line overrides the config file, check here if it is provided
if not force:
    force = args.force

force_assert = """
Stop, prefilter file already exists!
Use --force to overwrite.
Also check 'perform: true' is set in the config.
(Note that 'perform: false' will still output a file lacking prefiltered features.)
"""
if output_file.exists():
    assert force, force_assert

# Create the directory
output_file.parent.mkdir(exist_ok=True, parents=True)

# Perform prefiltering and output file
if perform:
    features_df = prefilter_features(core_args, example_site_dir, flag_cols)
else:
    features_df = load_features(core_args, example_dir)
    features_df = features_df.assign(prefilter_column=False)

features_df.to_csv(output_file, sep="\t", index=False)
