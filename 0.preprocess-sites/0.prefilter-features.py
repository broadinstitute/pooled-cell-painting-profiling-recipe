"""
0.prefilter-features.py

Determine which features should be used for building morphology profiles.

Note this is a preselection step.
An additional round of feature selection will occur at a different stage.
"""

import os
import sys
import pathlib
import warnings
import argparse
import numpy as np
import pandas as pd
from scripts.site_processing_utils import prefilter_features

sys.path.append(os.path.join("..", "scripts"))
from config_utils import process_config_file
from arg_utils import parse_command_args
from io_utils import check_if_write

args = parse_command_args(config_file="site_processing_config.yaml")
config_file = args.config_file
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

file_exist_warning = """
Warning, prefilter file already exists! Not overwriting!
Set 'force_overwrite: true' in config or use --force to overwrite.
Also check 'perform: true' is set in the config.
(Note that 'perform: false' will still output a file lacking prefiltered features.)
"""

force_warning = """
Warning, prefilter file already exists! Overwriting file. This may be intended.
"""

if output_file.exists():
    if not force:
        warnings.warn(file_exist_warning)
    else:
        warnings.warn(force_warning)

# Create the directory
output_file.parent.mkdir(exist_ok=True, parents=True)

# Perform prefiltering and output file
if perform:
    features_df = prefilter_features(core_args, example_site_dir, flag_cols)
else:
    features_df = load_features(core_args, example_dir)
    features_df = features_df.assign(prefilter_column=False)

if check_if_write(output_file, force):
    features_df.to_csv(output_file, sep="\t", index=False)
