"""
Does not respect overwrite
"""

import os
import sys
import json
import logging
import traceback
import pandas as pd

from pycytominer import aggregate

recipe_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(recipe_path, "utils"))
from io_utils import load_configs, read_csvs_with_chunksize, parse_data_set

# Configure logging
logfolder = os.path.join(recipe_path, "logs")
if not os.path.isdir(logfolder):
    os.mkdir(logfolder)
logging.basicConfig(
    filename=os.path.join(logfolder, "4.aggregate.log"),
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


def aggregate(path_to_defaults_config, path_to_experiment_config):
    printandlog("Starting 4.aggregate.")
    defaults_config, experiment_config = load_configs(
        path_to_defaults_config, path_to_experiment_config
    )

    # define experiment variables
    data_sets = experiment_config["data_sets"]

    # define defaults variables
    out_root = defaults_config["directory_structure"]["root"]
    single_cell = defaults_config["directory_structure"]["single_cell"]
    profiles = defaults_config["directory_structure"]["profiles"]
    data_folder = defaults_config["directory_structure"]["data"]
    compression = defaults_config["core"]["compression"]
    float_format = defaults_config["core"]["float_format"]
    barcode_col = defaults_config["process"]["process_SBS"]["barcode_col"]
    gene_col = defaults_config["process"]["process_SBS"]["gene_col"]
    save_single_file = defaults_config["process"]["single_cell"]["save_single_file"]
    allowed_skips = defaults_config["process"]["aggregate"]["allowed_skips"]
    operation = defaults_config["process"]["aggregate"]["operation"]
    features = defaults_config["process"]["aggregate"]["features"]

    outdir_data = os.path.join(out_root, data_folder)
    profiles_out = os.path.join(out_root, profiles)
    if not os.path.exists(profiles_out):
        os.makedirs(profiles_out, exist_ok=True)

    allowed_skip_counter = 0
    for data_set_name in data_sets.keys():
        printandlog(f"Starting {data_set_name}")
        batch, plate_list, list_plate_well_tuples = parse_data_set(
            data_sets[data_set_name]
        )
        folders_infile = os.path.join(outdir_data, f"{data_set_name}_sitelist.json")
        with open(folders_infile, "r") as file:
            usefolders = json.load(file)

        for plate in plate_list:
            # Load in single cell data
            if not save_single_file:
                single_file_df = []
                if allowed_skips >= allowed_skip_counter:
                    for plate_well_site_folder in [x for x in usefolders if plate in x]:
                        try:
                            infile = os.path.join(
                                out_root,
                                single_cell,
                                "by_site",
                                plate_well_site_folder,
                                f"{plate_well_site_folder}_single_cell.csv.gz",
                            )
                            df = read_csvs_with_chunksize(infile)
                            single_file_df.append(df)
                        except:
                            printandlog(
                                f"Skipped loading per-site single cell data for {plate_well_site_folder}.",
                                type="warning",
                            )
                            allowed_skip_counter += 1
                            printandlog(
                                f"Now at {allowed_skip_counter} sites skipped from errors.",
                                type="warning",
                            )
                single_file_df = pd.concat(single_file_df, axis="rows").reset_index(
                    drop=True
                )
            else:
                infile = os.path.join(
                    out_root, single_cell, f"{plate}_single_cell_{data_set_name}.csv.gz"
                )
                single_file_df = read_csvs_with_chunksize(infile)

            # Aggregate to guide level
            printandlog(
                f"Now aggregating {plate} to guide level with operation: {operation}"
            )

            aggregate_df = aggregate(
                population_df=single_file_df,
                strata=[f"Metadata_{gene_col}", f"Metadata_{barcode_col}"],
                features=features,
                operation=operation,
            )
            filepath = os.path.join(
                profiles_out, f"{plate}_{data_set_name}_guide.csv.gz"
            )
            aggregate_df.to_csv(
                filepath,
                compression=compression,
                float_format=float_format,
                index=False,
            )

            # Aggregate to gene level
            printandlog(f"Now aggregating to gene level with operation: {operation}")

            aggregate_df = aggregate(
                population_df=single_file_df,
                strata=[f"Metadata_{gene_col}"],
                features=features,
                operation=operation,
            )
            filepath = os.path.join(
                profiles_out, f"{plate}_{data_set_name}_gene.csv.gz"
            )
            aggregate_df.to_csv(
                filepath,
                compression=compression,
                float_format=float_format,
                index=False,
            )
    printandlog("Finished 4.aggregate.")


#################################
# MAIN USER INTERACTION
#################################

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Pass paths to defaults config and experiment config.")
        sys.exit()
    path_to_defaults_config = sys.argv[1]
    path_to_experiment_config = sys.argv[2]
    aggregate(path_to_defaults_config, path_to_experiment_config)
