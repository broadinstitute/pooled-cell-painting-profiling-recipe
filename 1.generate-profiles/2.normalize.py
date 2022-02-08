import os
import sys
import pathlib
import argparse
import warnings
import logging
import traceback
import pandas as pd
from joblib import Parallel, delayed

from pycytominer import normalize

sys.path.append("config")
from utils import parse_command_args, process_configuration, get_split_aware_site_info

recipe_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(recipe_path, "scripts"))
from io_utils import read_csvs_with_chunksize

# Configure logging
logfolder = os.path.join(os.path.dirname(recipe_path), "logs")
if not os.path.isdir(logfolder):
    os.mkdir(logfolder)
logging.basicConfig(
    filename=os.path.join(logfolder, "2.normalize.log"), level=logging.INFO,
)


def handle_excepthook(exc_type, exc_value, exc_traceback):
    logging.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))
    traceback_details = "\n".join(traceback.extract_tb(exc_traceback).format())
    print(f"Uncaught Exception: {traceback_details}")


sys.excepthook = handle_excepthook

# Configure experiment
args = parse_command_args()
logging.info(f"Args used:{args}")

batch_id = args.batch_id
options_config_file = args.options_config_file
experiment_config_file = args.experiment_config_file
split_step = args.split_step

config, incomplete_sites, errored_sites = process_configuration(
    batch_id,
    step="profile--normalize",
    options_config=options_config_file,
    experiment_config=experiment_config_file,
)
logging.info(f"Config used:{config}")
logging.info(f"Skipped incomplete sites during config processing: {incomplete_sites}")
logging.info(f"Skipped errored sites during config processing: {errored_sites}")

# Extract config arguments
split_info = config["experiment"]["split"][split_step]
perform = config["options"]["profile"]["normalize"]["perform"]

# check if this step should be performed
if not perform:
    sys.exit("Config file set to perform=False, not performing {}".format(__file__))

ignore_files = config["options"]["core"]["ignore_files"]
float_format = config["options"]["core"]["float_format"]
compression = config["options"]["core"]["compression"]

input_spotdir = config["directories"]["preprocess"]["spots"]
normalize_input_dir = config["directories"]["profile"]["profiles"]
single_cell_input_dir = config["directories"]["profile"]["single_cell"]
normalize_input_files = config["files"]["aggregate_files"]
normalize_output_files = config["files"]["normalize_files"]
single_cell_file = config["files"]["single_file_only_output_file"]
image_file = config["files"]["image_file"]

sc_config = config["options"]["profile"]["single_cell"]
normalize_singlecell_from_single_file = sc_config["output_one_single_cell_file_only"]

normalize_args = config["options"]["profile"]["normalize"]
output_single_cell_by_guide = normalize_args["output_single_cell_by_guide"]
normalize_levels = normalize_args["levels"]
normalize_by_samples = normalize_args["by_samples"]
normalize_these_features = normalize_args["features"]
normalize_method = normalize_args["method"]
force = normalize_args["force_overwrite"]

print("Starting 2.normalize.")
logging.info(f"Started 2.normalize.")

sites = [x.name for x in input_spotdir.iterdir() if x.name not in ignore_files]
site_info_dict = get_split_aware_site_info(
    config["experiment"], sites, split_info, separator="___"
)

def append_to_guide_csv(guide, image_df):
    append_df = df.loc[
        df["Metadata_Foci_Barcode_MatchedTo_Barcode"] == guide
    ]
    gene = list(append_df["Metadata_Foci_Barcode_MatchedTo_Barcode"])[0]
    guide_file_name = f"{str(output_file).split('__')[0].split('/')[-1]}__{guide}_{gene}.csv.gz"
    guide_path = os.path.join(sc_by_guide_folder, guide_file_name)
    if not os.path.exists(guide_path):
        guide_df = pd.DataFrame()
    else:
        guide_df = read_csvs_with_chunksize(guide_path)

    append_df = append_df.merge(
        image_df, left_on="Metadata_Foci_site", right_on="Metadata_site"
    )
    guide_df = guide_df.append(append_df)
    guide_df.to_csv(guide_path, index=False)


for data_split_site in site_info_dict:
    for data_level in normalize_levels:
        if data_level == "single_cell":
            if not normalize_singlecell_from_single_file:
                continue
            file_to_normalize = pathlib.Path(
                single_cell_input_dir,
                single_cell_file.name.replace(".csv.gz", f"_{data_split_site}.csv.gz"),
            )
        else:
            file_to_normalize = normalize_input_files[data_level]
            file_to_normalize = pathlib.Path(
                normalize_input_dir,
                file_to_normalize.name.replace(".csv.gz", f"_{data_split_site}.csv.gz"),
            )

        print(
            f"Now normalizing {data_level}...with operation: {normalize_method} for split {data_split_site}"
        )
        logging.info(
            f"Normalizing {data_level}...with operation: {normalize_method} for split {data_split_site}"
        )

        output_file = normalize_output_files[data_level]
        output_file = pathlib.Path(
            normalize_output_files[data_level].parents[0],
            output_file.name.replace(".csv.gz", f"_{data_split_site}.csv.gz"),
        )
        df = read_csvs_with_chunksize(file_to_normalize)

        normalize(
            profiles=df,
            features=normalize_these_features,
            samples=normalize_by_samples,
            method=normalize_method,
            output_file=output_file,
            compression_options=compression,
            float_format=float_format,
        )

        if data_level == "single_cell":
            if output_single_cell_by_guide:
                print(
                    f"Now outputting normalized single cell profiles by guide for split {data_split_site}"
                )
                logging.info(
                    f"Now outputting normalized single cell profiles by guide for split {data_split_site}"
                )
                # Load image alignment information for appending to single_cell_by_guide csvs
                image_df = pd.read_csv(image_file, sep="\t")
                keep_columns = []
                for col in image_df.columns:
                    if "Align_" in col:
                        keep_columns.append(col)
                keep_columns.append("Metadata_site")
                image_df = image_df.loc[:, keep_columns]

                sc_by_guide_folder = os.path.join(
                    single_cell_input_dir, "single_cell_by_guide"
                )
                if not os.path.isdir(sc_by_guide_folder):
                    os.mkdir(sc_by_guide_folder)
                df = read_csvs_with_chunksize(output_file)
                Parallel(n_jobs=-2)(delayed(append_to_guide_csv)(guide, image_df) for guide in set(df["Metadata_Foci_Barcode_MatchedTo_Barcode"]))

print("Finished 2.normalize.")
logging.info("Finished 2.normalize.")
