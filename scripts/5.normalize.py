"""
Outputs the following files:
- {profiles}/{plate}_{data_set_name}_{gene|guide}_normalized.csv.gz
- {single_cell}/single_cell_by_gene/{data_set_name}_gene_normalized_{gene}.csv.gz  (single cells, only if output_bygene)
- {single_cell}/single_cell_by_gene/{data_set_name}_guide_normalized_{barcode}.csv.gz  (single cells, only if output_bygene)
"""

import os
import sys
import glob
import json
import logging
import traceback
import pandas as pd

from pycytominer import normalize as pycytominer_normalize

recipe_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(recipe_path, "utils"))
from io_utils import (
    load_configs,
    read_csvs_with_chunksize,
    parse_data_set,
    site_folders_for_plate,
)

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


def prepare_for_normalization(df, by_samples, control_genes, gene_col):
    """
    Work out which columns are metadata/locations (not normalized) and add the
    Metadata_Norm column that selects the samples normalization is calculated from.
    Returns (df, meta_cols), or (None, None) if by_samples couldn't be parsed.
    """
    # Don't normalize locations
    meta_cols = list(df.columns[df.columns.str.contains("Metadata")])
    remove_locs = list(
        filter(
            lambda x: "_Location_Center_X" in x or "_Location_Center_Y" in x,
            df.columns,
        )
    )
    remove_cents = list(
        filter(
            lambda x: "AreaShape_Center_X" in x or "AreaShape_Center_Y" in x,
            df.columns,
        )
    )
    meta_cols = meta_cols + remove_locs + remove_cents

    # Set which samples to normalize to by making Metadata_Norm column
    if by_samples == "control_genes":
        df["Metadata_Norm"] = "False"
        df.loc[df[f"Metadata_{gene_col}"].isin(control_genes), "Metadata_Norm"] = "True"
    elif by_samples == "all":
        df["Metadata_Norm"] = "True"
    else:
        return None, None
    return df, meta_cols


def load_single_cell_plate(
    out_root,
    single_cell_folder,
    data_set_name,
    plate,
    save_single_file,
    usefolders,
):
    """Load all single cells of a plate, from the per-site files or the single plate file."""
    if save_single_file:
        return read_csvs_with_chunksize(
            os.path.join(
                out_root,
                single_cell_folder,
                f"{plate}_single_cell_{data_set_name}.csv.gz",
            )
        )
    dfs = []
    for plate_well_site_folder in site_folders_for_plate(usefolders, plate):
        infile = os.path.join(
            out_root,
            single_cell_folder,
            "by_site",
            plate_well_site_folder,
            f"{plate_well_site_folder}_single_cell.csv.gz",
        )
        try:
            site_df = read_csvs_with_chunksize(infile)
            # sites without cells would turn numeric columns into object dtype
            if len(site_df) > 0:
                dfs.append(site_df)
        except Exception:
            printandlog(
                f"Skipped {plate_well_site_folder}. Couldn't load {infile}",
                type="warning",
            )
    if not dfs:
        return None
    return pd.concat(dfs, axis="rows").reset_index(drop=True)


def append_by_group(df, group_col, name, sc_by_gene_folder, data_set_name, compression, float_format):
    """Write each group's single cells to {data_set}_{name}_normalized_{group}.csv.gz, appending across plates."""
    for group, slice in df.groupby(group_col):
        group_path = os.path.join(
            sc_by_gene_folder, f"{data_set_name}_{name}_normalized_{group}.csv.gz"
        )
        if os.path.exists(group_path):
            slice = pd.concat([read_csvs_with_chunksize(group_path), slice])
        slice.to_csv(
            group_path,
            index=False,
            compression=compression,
            float_format=float_format,
        )


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
    data_folder = defaults_config["directory_structure"]["data"]
    gene_col = defaults_config["process"]["process_SBS"]["gene_col"]
    barcode_col = defaults_config["process"]["process_SBS"]["barcode_col"]
    save_single_file = defaults_config["process"]["single_cell"]["save_single_file"]
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

        sc_by_gene_folder = os.path.join(
            out_root, single_cell_folder, "single_cell_by_gene"
        )
        if output_bygene:
            # by-gene files are appended to plate by plate, so clear previous runs' files first
            for normby in ["gene", "guide"]:
                for old_file in glob.glob(
                    os.path.join(
                        glob.escape(sc_by_gene_folder),
                        f"{glob.escape(data_set_name)}_{normby}_normalized_*.csv.gz",
                    )
                ):
                    os.remove(old_file)
            usefolders = []
            if not save_single_file:
                with open(
                    os.path.join(
                        out_root, data_folder, f"{data_set_name}_sitelist.json"
                    ),
                    "r",
                ) as file:
                    usefolders = json.load(file)

        for plate in plate_list:
            for normby in ["gene", "guide"]:
                printandlog(f"Normalizing {data_set_name} {plate} by {normby}.")
                df = read_csvs_with_chunksize(
                    os.path.join(
                        profiles_out, f"{plate}_{data_set_name}_{normby}.csv.gz"
                    )
                )
                df = df.loc[:, (df != 0).any(axis=0)]  # drop columns with all 0

                df, meta_cols = prepare_for_normalization(
                    df, by_samples, control_genes, gene_col
                )
                if df is None:
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

            if output_bygene:
                printandlog(
                    f"Normalizing single cells of {data_set_name} {plate} for by-gene output."
                )
                df = load_single_cell_plate(
                    out_root,
                    single_cell_folder,
                    data_set_name,
                    plate,
                    save_single_file,
                    usefolders,
                )
                if df is None:
                    printandlog(
                        f"No single cell data loaded for {data_set_name} {plate}. Skipping by-gene output.",
                        type="warning",
                    )
                    continue
                # drop all-0 feature columns, but never metadata
                df = df.loc[
                    :,
                    [c for c in df.columns if "Metadata" in c or (df[c] != 0).any()],
                ]
                df, meta_cols = prepare_for_normalization(
                    df, by_samples, control_genes, gene_col
                )
                if df is None:
                    printandlog(
                        "Failed to parse appropriate normalization by_samples. Check your by_samples value.",
                        type="warning",
                    )
                    return
                if len(df.loc[df["Metadata_Norm"] == "True"]) == 0:
                    printandlog(
                        f"No samples were selected for normalization for {data_set_name} {plate} single cells. Check your control_genes and that your sample contains the control.",
                        type="warning",
                    )
                    continue
                df = pycytominer_normalize(
                    profiles=df,
                    features=features,
                    meta_features=meta_cols,
                    samples="Metadata_Norm == 'True'",
                    method=method,
                )
                os.makedirs(sc_by_gene_folder, exist_ok=True)
                append_by_group(
                    df, f"Metadata_{gene_col}", "gene", sc_by_gene_folder,
                    data_set_name, compression, float_format,
                )
                append_by_group(
                    df, f"Metadata_{barcode_col}", "guide", sc_by_gene_folder,
                    data_set_name, compression, float_format,
                )

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
