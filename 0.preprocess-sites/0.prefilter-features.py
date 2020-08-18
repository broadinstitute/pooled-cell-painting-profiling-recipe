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

sys.path.append("config")
from utils import parse_command_args, process_configuration

recipe_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(recipe_path, "scripts"))
from io_utils import check_if_write

args = parse_command_args()

plate_id = args.plate_id
options_config_file = args.options_config_file
experiment_config_file = args.experiment_config_file

config = process_configuration(
    plate_id,
    options_config=options_config_file,
    experiment_config=experiment_config_file,
)

# Set constants
experiment_args = config["experiment"]
prefilter_args = config["options"]["preprocess"]["prefilter"]
core_option_args = config["options"]["core"]
force = prefilter_args["force_overwrite"]
perform = prefilter_args["perform"]
flag_cols = prefilter_args["flag_cols"]

prefilter_file = config["files"]["prefilter_file"]
input_dir = config["directories"]["input_data_dir"]
example_site = config["options"]["example_site"]
example_site_dir = pathlib.Path(f"{input_dir}/{example_site}")

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

if prefilter_file.exists():
    if not force:
        warnings.warn(file_exist_warning)
    else:
        warnings.warn(force_warning)

# Perform prefiltering and output file
if perform:
    features_df = prefilter_features(core_option_args, example_site_dir, flag_cols)
else:
    features_df = load_features(core_option_args, example_site_dir)
    features_df = features_df.assign(prefilter_column=False)

if check_if_write(prefilter_file, force):
    features_df.to_csv(prefilter_file, sep="\t", index=False)
