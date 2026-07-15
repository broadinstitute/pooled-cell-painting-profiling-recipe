"""
Outputs the following files:
- {profiles}/{plate}_{data_set_name}_{gene|guide}_normalized.csv.gz
- {single_cell}/{single_cell}/single_cell_by_gene/{data_set_name}_{gene|guide}_normalized_{gene}_{gene}.csv.gz
"""

import os
import sys
import logging
import traceback
import pandas as pd

from pycytominer import normalize as pycytominer_normalize

recipe_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(recipe_path, "utils"))
from io_utils import load_configs, read_csvs_with_chunksize, parse_data_set

# Configure logging
logfolder = os.path.join(recipe_path, "logs")
if not os.path.isdir(logfolder):
    os.mkdir(logfolder)
logging.basicConfig(
    filename=os.path.join(logfolder, "5.normalize.log"),
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


def normalize(path_to_defaults_config, path_to_experiment_config):
    printandlog("Starting 5.normalize.")
    defaults_config, experiment_config = load_configs(
        path_to_defaults_config, path_to_experiment_config
    )

    # define experiment variables
    data_sets = experiment_config["data_sets"]
    control_genes = experiment_config["control_genes"]

    # define defaults variables
    out_root = defaults_config["directory_structure"]["root"]
    profiles = defaults_config["directory_structure"]["profiles"]
    single_cell_folder = defaults_config["directory_structure"]["single_cell"]
    compression = defaults_config["core"]["compression"]
    float_format = defaults_config["core"]["float_format"]
    gene_col = defaults_config["process"]["process_SBS"]["gene_col"]
    output_bygene = defaults_config["process"]["normalize"]["output_bygene"]
    method = defaults_config["process"]["normalize"]["method"]
    by_samples = defaults_config["process"]["normalize"]["by_samples"]
    features = defaults_config["process"]["normalize"]["features"]

    profiles_out = os.path.join(out_root, profiles)

    for data_set_name in data_sets.keys():
        printandlog(f"Starting {data_set_name}")
        batch, plate_list, list_plate_well_tuples = parse_data_set(
            data_sets[data_set_name]
        )

        for plate in plate_list:
            for normby in ["gene", "guide"]:
                skip_group = False
                printandlog(f"Normalizing {data_set_name} {plate} by {normby}.")
                df = read_csvs_with_chunksize(
                    os.path.join(
                        profiles_out, f"{plate}_{data_set_name}_{normby}.csv.gz"
                    )
                )
                df = df.loc[:, (df != 0).any(axis=0)]  # drop columns with all 0

                # Don't normalize locations
                meta_cols = list(df.columns[df.columns.str.contains("Metadata")])
                remove_locs = list(
                    filter(
                        lambda x: "_Location_Center_X" in x
                        or "_Location_Center_Y" in x,
                        df.columns,
                    )
                )
                remove_cents = list(
                    filter(
                        lambda x: "AreaShape_Center_X" in x
                        or "AreaShape_Center_Y" in x,
                        df.columns,
                    )
                )
                meta_cols = meta_cols + remove_locs + remove_cents

                # Set which samples to normalize to by making Metadata_Norm column
                if by_samples == "control_genes":
                    df["Metadata_Norm"] = "False"
                    df.loc[
                        df[f"Metadata_{gene_col}"].isin(control_genes), "Metadata_Norm"
                    ] = "True"
                elif by_samples == "all":
                    df["Metadata_Norm"] = "True"
                else:
                    printandlog(
                        "Failed to parse appropriate normalization by_samples. Check your by_samples value.",
                        type="warning",
                    )
                    return
                if len(df.loc[df["Metadata_Norm"] == "True"]) == 0:
                    printandlog(
                        f"No samples were selected for normalization for {data_set_name} {plate} {normby}. Check your control_genes and that your sample contains the control.",
                        type="warning",
                    )
                    skip_group = True
                    continue
                output_file = os.path.join(
                    profiles_out, f"{plate}_{data_set_name}_{normby}_normalized.csv.gz"
                )
                pycytominer_normalize(
                    profiles=df,
                    features=features,
                    meta_features=meta_cols,
                    samples="Metadata_Norm == 'True'",
                    method=method,
                    output_file=output_file,
                    compression_options=compression,
                    float_format=float_format,
                )
            if output_bygene and not skip_group:
                sc_by_gene_folder = os.path.join(
                    out_root, single_cell_folder, "single_cell_by_gene"
                )
                if not os.path.isdir(sc_by_gene_folder):
                    os.mkdir(sc_by_gene_folder)

                df = read_csvs_with_chunksize(output_file)
                for gene in df[f"Metadata_{gene_col}"].unique():
                    slice = df.loc[df[f"Metadata_{gene_col}"] == gene]
                    gene_file_name = (
                        f"{data_set_name}_{normby}_normalized_{gene}_{gene}.csv.gz"
                    )
                    gene_path = os.path.join(sc_by_gene_folder, gene_file_name)
                    if not os.path.exists(gene_path):
                        gene_df = pd.DataFrame()
                    else:
                        gene_df = read_csvs_with_chunksize(gene_path)
                    gene_df = pd.concat([gene_df, slice])
                    gene_df.to_csv(gene_path, index=False)

    printandlog("Finished 5.normalize.")


#################################
# MAIN USER INTERACTION
#################################

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Pass paths to defaults config and experiment config.")
        sys.exit()
    path_to_defaults_config = sys.argv[1]
    path_to_experiment_config = sys.argv[2]
    normalize(path_to_defaults_config, path_to_experiment_config)
