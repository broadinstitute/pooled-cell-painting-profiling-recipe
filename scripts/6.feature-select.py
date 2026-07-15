"""
Outputs the following files:
- profiles/{plate}_{data_set_name}_{gene|guide}_normalized_{group}_feature_selected.csv.gz
"""

import os
import sys
import warnings
import logging
import traceback
import pandas as pd

from pycytominer import feature_select as pycytominer_feature_select

recipe_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(recipe_path, "utils"))
from io_utils import load_configs, read_csvs_with_chunksize, parse_data_set

# Configure logging
logfolder = os.path.join(recipe_path, "logs")
if not os.path.isdir(logfolder):
    os.mkdir(logfolder)
logging.basicConfig(
    filename=os.path.join(logfolder, "6.feature-select.log"),
    level=logging.INFO,
)


def handle_excepthook(exc_type, exc_value, exc_traceback):
    logging.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))
    traceback_details = "\n".join(traceback.extract_tb(exc_traceback).format())
    print(f"Uncaught Exception: {traceback_details}")


sys.excepthook = handle_excepthook


def printandlog(msg, type="info"):
    print(msg)
    if type == "warning":
        logging.warning(msg)
    else:
        logging.info(msg)


def feature_select(path_to_defaults_config, path_to_experiment_config):
    printandlog("Starting 6.feature_select.")
    defaults_config, experiment_config = load_configs(
        path_to_defaults_config, path_to_experiment_config
    )

    # define experiment variables
    data_sets = experiment_config["data_sets"]

    # define defaults variables
    out_root = defaults_config["directory_structure"]["root"]
    profiles = defaults_config["directory_structure"]["profiles"]
    compression = defaults_config["core"]["compression"]
    float_format = defaults_config["core"]["float_format"]
    group = defaults_config["process"]["feature_select"]["group"]
    operations = defaults_config["process"]["feature_select"]["operations"]
    use_samples = defaults_config["process"]["feature_select"]["use_samples"]
    features = defaults_config["process"]["feature_select"]["features"]
    na_cutoff = defaults_config["process"]["feature_select"]["na_cutoff"]
    corr_threshold = defaults_config["process"]["feature_select"]["corr_threshold"]

    profiles_out = os.path.join(out_root, profiles)

    for normby in ["gene", "guide"]:
        all_normed = pd.DataFrame()
        for data_set_name in data_sets.keys():
            printandlog(f"Starting {data_set_name} {normby}")
            batch, plate_list, list_plate_well_tuples = parse_data_set(
                data_sets[data_set_name]
            )
            group_normed = pd.DataFrame()
            for plate in plate_list:
                df = read_csvs_with_chunksize(
                    os.path.join(
                        profiles_out,
                        f"{plate}_{data_set_name}_{normby}_normalized.csv.gz",
                    )
                )
                if group == "plate":
                    pycytominer_feature_select(
                        profiles=df,
                        features=features,
                        samples=use_samples,
                        operation=operations,
                        na_cutoff=na_cutoff,
                        corr_threshold=corr_threshold,
                        output_file=os.path.join(
                            out_root,
                            f"{plate}_{data_set_name}_{normby}_normalized_{group}_feature_selected.csv.gz",
                        ),
                        compression_options=compression,
                        float_format=float_format,
                    )
                else:
                    group_normed = pd.concat([group_normed, df])
            if group == "group":
                pycytominer_feature_select(
                    profiles=group_normed,
                    features=features,
                    samples=use_samples,
                    operation=operations,
                    na_cutoff=na_cutoff,
                    corr_threshold=corr_threshold,
                    output_file=os.path.join(
                        out_root,
                        f"{plate}_{data_set_name}_{normby}_normalized_{group}_feature_selected.csv.gz",
                    ),
                    compression_options=compression,
                    float_format=float_format,
                )
            elif group == "all":
                all_normed = pd.concat([all_normed, df])
        pycytominer_feature_select(
            profiles=all_normed,
            features=features,
            samples=use_samples,
            operation=operations,
            na_cutoff=na_cutoff,
            corr_threshold=corr_threshold,
            output_file=os.path.join(
                out_root,
                f"{plate}_{data_set_name}_{normby}_normalized_{group}_feature_selected.csv.gz",
            ),
            compression_options=compression,
            float_format=float_format,
        )
    printandlog("Done with 6.feature_select")


#################################
# MAIN USER INTERACTION
#################################

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Pass paths to defaults config and experiment config.")
        sys.exit()
    path_to_defaults_config = sys.argv[1]
    path_to_experiment_config = sys.argv[2]
    feature_select(path_to_defaults_config, path_to_experiment_config)
