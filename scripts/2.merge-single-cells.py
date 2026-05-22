"""
Allows for resuming and either overwriting or not.
List of folders used is saved out for future steps.
Last step that files are read from CellProfiler outputs.
"""

import os
import sys
import warnings
import logging
import json
import traceback
import pandas as pd

recipe_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(recipe_path, "utils"))
from io_utils import load_configs, parse_data_set, read_csvs_with_chunksize
from cell_quality_utils import get_cell_quality_dict, filter_to_top_BC

# Configure logging
logfolder = os.path.join(recipe_path, "logs")
if not os.path.isdir(logfolder):
    os.mkdir(logfolder)
logging.basicConfig(
    filename=os.path.join(logfolder, "2.merge-single-cells.log"),
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


def fail_site(plate_well_site_folder, file, allowed_skip_counter, allowed_skips):
    printandlog(
        f"Skipped {plate_well_site_folder}. Couldn't find or load necessary columns from {file}",
        type="warning",
    )
    allowed_skip_counter += 1
    printandlog(
        f"Now at {allowed_skip_counter} sites skipped from errors.",
        type="warning",
    )
    if allowed_skips <= allowed_skip_counter:
        printandlog(
            f"Allowed skip limit of {allowed_skips} reached. Stopping 2.merge-single-cells.",
            type="warning",
        )
        sys.exit(1)
    else:
        return allowed_skip_counter


def merge_single_cells(path_to_defaults_config, path_to_experiment_config):
    printandlog("Starting 2.merge-single-cells.")
    defaults_config, experiment_config = load_configs(
        path_to_defaults_config, path_to_experiment_config
    )

    # define experiment variables
    data_sets = experiment_config["data_sets"]
    file_location = experiment_config["file_location"]
    compartments = experiment_config["compartments"]
    cell_quality_method = experiment_config["cell_quality_method"]
    overwrite_files = experiment_config["overwrite_files"]

    # define defaults variables
    out_root = defaults_config["directory_structure"]["root"]
    sbs_folder = defaults_config["directory_structure"]["SBS"]
    data_folder = defaults_config["directory_structure"]["data"]
    single_cell = defaults_config["directory_structure"]["single_cell"]
    cell_quality_column = defaults_config["core"]["cell_quality_column"]
    cell_quality_index = defaults_config["core"]["cell_quality_index"]
    cell_id_cols = defaults_config["core"]["cell_id_cols"]
    compression = defaults_config["core"]["compression"]
    float_format = defaults_config["core"]["float_format"]
    SBS_score_col = defaults_config["process"]["process_SBS"]["SBS_score_col"]
    allowed_skips = defaults_config["process"]["single_cell"]["allowed_skips"]
    flag_cols = defaults_config["process"]["single_cell"]["flag_cols"]
    save_single_file = defaults_config["process"]["single_cell"]["save_single_file"]
    save_per_site_file = defaults_config["process"]["single_cell"]["save_per_site_file"]
    cell_quality_dict = get_cell_quality_dict(cell_quality_method)

    single_cell_out = os.path.join(out_root, single_cell, "by_site")
    if not os.path.exists(single_cell_out):
        os.makedirs(single_cell_out, exist_ok=True)
    outdir_data = os.path.join(out_root, data_folder)
    if not os.path.exists(outdir_data):
        os.makedirs(outdir_data, exist_ok=True)

    if not overwrite_files:
        printandlog(f"Will not overwrite existing single-cell files.")
    if save_single_file and save_per_site_file:
        printandlog(f"Will save out a file per site AND a single file for all sites.")
        if not overwrite_files:
            printandlog(
                f"This may lead to unexpected behavior in single-file because it will NOT include previously generated sites"
            )
    elif save_single_file:
        printandlog(
            f"Will save out only a single file for all sites and NOT a file per site."
        )
    elif save_per_site_file:
        printandlog(
            f"Will save out only a file per site and NOT a single file for all sites."
        )

    allowed_skip_counter = 0
    for data_set_name in data_sets.keys():
        printandlog(f"Starting {data_set_name}")
        batch, plate_list, list_plate_well_tuples = parse_data_set(
            data_sets[data_set_name]
        )

        folderlist = os.listdir(os.path.join(file_location, batch))
        try:
            inferred_empty_sites = (
                pd.read_csv(
                    os.path.join(
                        out_root,
                        data_folder,
                        f"Inferred_Empty_Sites_{data_set_name}.csv",
                    )
                )
                .squeeze()
                .tolist()
            )
        except:
            inferred_empty_sites = []
            printandlog(f"No inferred empty sites file found for {data_set_name}.")
        try:
            no_SBS_foci_sites = (
                pd.read_csv(
                    os.path.join(
                        out_root,
                        data_folder,
                        f"No_SBS_Foci_Sites_{data_set_name}.csv",
                    )
                )
                .squeeze()
                .tolist()
            )
        except:
            no_SBS_foci_sites = []
            printandlog(f"No SBS foci-less file found for {data_set_name}.")

        all_usefolders = []
        for plate in plate_list:
            for plate_well_tuple in [x for x in list_plate_well_tuples if plate in x]:
                well = plate_well_tuple[1]
                usefolders = [x for x in folderlist if f"{plate}-{well}" in x]
                usefolders = [x for x in usefolders if x not in inferred_empty_sites]
                usefolders = [x for x in usefolders if x not in no_SBS_foci_sites]
                all_usefolders += usefolders
                single_file_df = []
                for plate_well_site_folder in usefolders:
                    plate_well_site_out = os.path.join(
                        single_cell_out, plate_well_site_folder
                    )
                    SBS_out = os.path.join(out_root, sbs_folder, plate_well_site_folder)
                    if save_per_site_file:
                        if not overwrite_files:
                            if os.path.exists(plate_well_site_out):
                                printandlog(
                                    f"Skipping {plate_well_site_folder}. Output folder {plate_well_site_out} exists and overwrite_files is False."
                                )
                                continue

                    site = plate_well_site_folder.rsplit("-", 1)[1]

                    # Load processed SBS data
                    indir = os.path.join(out_root, sbs_folder, plate_well_site_folder)
                    try:
                        SBS_df = read_csvs_with_chunksize(
                            os.path.join(indir, "processed_SBS.tsv.gz"),
                            sep="\t",
                        )
                    except:
                        allowed_skip_counter = fail_site(
                            plate_well_site_folder,
                            "processed_SBS.tsv.gz",
                            allowed_skip_counter,
                            allowed_skips,
                        )
                        continue

                    # Filter SBS to top quality barcode/Cell
                    SBS_df = filter_to_top_BC(SBS_df, compartments[0], [x for x in SBS_df.columns if SBS_score_col in x and '_mean' in x][0])
                    # Rename columns in preparation for merge
                    SBS_df.columns = [
                        f"Metadata_{x}" if "Metadata" not in x else x
                        for x in SBS_df.columns
                    ]

                    # Load csv files for prespecified compartments
                    indir = os.path.join(file_location, batch, plate_well_site_folder)
                    compartment_csvs = {}
                    try:
                        for compartment in compartments:
                            compartment_df = read_csvs_with_chunksize(
                                os.path.join(indir, f"{compartment}.csv")
                            )
                            # Remove any columns that have any string from flag_cols in them
                            compartment_df = compartment_df[
                                [
                                    x
                                    for x in compartment_df.columns
                                    if not any(flag in x for flag in flag_cols)
                                ]
                            ]
                            # Rename columns in preparation for merging
                            compartment_df.columns = [
                                f"{compartment}_{x}" if x not in cell_id_cols else x
                                for x in compartment_df.columns
                            ]
                            compartment_df.columns = [
                                f"Metadata_{x}" if x in cell_id_cols else x
                                for x in compartment_df.columns
                            ]
                            compartment_csvs[compartment] = compartment_df
                    except FileNotFoundError:
                        allowed_skip_counter = fail_site(
                            plate_well_site_folder,
                            f"{compartment}.csv",
                            allowed_skip_counter,
                            allowed_skips,
                        )
                        continue

                    # Start with the first DataFrame
                    sc_merged_df = compartment_csvs[compartments[0]]
                    # Merge with subsequent DataFrames
                    for key, df in compartment_csvs.items():
                        if key != compartments[0]:  # Skip the first DataFrame
                            sc_merged_df = pd.merge(
                                sc_merged_df,
                                df,
                                on=[f"Metadata_{x}" for x in cell_id_cols],
                                how="outer",
                            )

                    assert (
                        len(sc_merged_df) > 0
                    ), "Merge failed to create a DataFrame with rows!"

                    sc_merged_df = sc_merged_df.assign(
                        Metadata_Batch=batch,
                        Metadata_Plate=plate,
                        Metadata_Well=well,
                        Metadata_Site=site,
                        Metadata_Identifier=f"{batch}-{plate}-{well}-{site}",
                        Metadata_Dataset_Split=data_set_name,
                    )

                    sc_merged_df.Metadata_Site = sc_merged_df.Metadata_Site.astype(int)
                    SBS_df.Metadata_Site = SBS_df.Metadata_Site.astype(int)

                    # Merge SBS data into phenotyping data
                    sc_merged_df = sc_merged_df.merge(
                        SBS_df,
                        left_on=[
                            "Metadata_ObjectNumber",
                            "Metadata_ImageNumber",
                            "Metadata_Batch",
                            "Metadata_Plate",
                            "Metadata_Well",
                            "Metadata_Site",
                            "Metadata_Identifier",
                            "Metadata_Dataset_Split",
                        ],
                        right_on=[
                            f"Metadata_Parent_{compartments[0]}",
                            "Metadata_ImageNumber",
                            "Metadata_Batch",
                            "Metadata_Plate",
                            "Metadata_Well",
                            "Metadata_Site",
                            "Metadata_Identifier",
                            "Metadata_Dataset_Split",
                        ],
                        how="left",
                    )

                    assert (
                        len(sc_merged_df) > 0
                    ), "Merge with SBS failed to create a DataFrame with rows!"

                    # Create new Quality of "Empty" for cells without SBS info
                    sc_merged_df.loc[:, cell_quality_column] = sc_merged_df.loc[
                        :, cell_quality_column
                    ].fillna("Empty")
                    sc_merged_df.loc[:, cell_quality_index] = sc_merged_df.loc[
                        :, cell_quality_index
                    ].fillna(len(cell_quality_dict) + 1)

                    if not os.path.exists(SBS_out):
                        os.makedirs(SBS_out, exist_ok=True)

                    cols = [
                        "Metadata_Plate",
                        "Metadata_Well",
                        "Metadata_Site",
                        "Metadata_Dataset_Split",
                        "Metadata_Quality_Name",
                    ]
                    plate_well_site_quality_summary = (
                        sc_merged_df[cols + ["Metadata_ObjectNumber"]]
                        .groupby(cols)
                        .count()
                        .reset_index()
                        .rename(columns={"Metadata_ObjectNumber": "count"})
                    )
                    plate_well_site_quality_summary.to_csv(
                        os.path.join(
                            SBS_out,
                            f"{plate_well_site_folder}_cell_quality_summary.csv",
                        ),
                        index=False,
                    )
                    if save_single_file:
                        single_file_df.append(sc_merged_df)
                    if save_per_site_file:
                        if not os.path.exists(plate_well_site_out):
                            os.makedirs(plate_well_site_out, exist_ok=True)
                        sc_merged_df.to_csv(
                            os.path.join(
                                plate_well_site_out,
                                f"{plate_well_site_folder}_single_cell.csv.gz",
                            ),
                            sep=",",
                            index=False,
                            compression=compression,
                            float_format=float_format,
                        )
        if save_single_file:
            outfile = os.path.join(
                out_root, single_cell, f"{plate}_single_cell_{data_set_name}.csv.gz"
            )
            if not overwrite_files:
                if os.path.exists(outfile):
                    printandlog(
                        f"Did not save single file for {data_set_name} {plate}. File already exists and overwrite files is False"
                    )
                    continue
            pd.concat(single_file_df, axis="rows").reset_index(drop=True).to_csv(
                outfile,
                sep=",",
                index=False,
                compression=compression,
                float_format=float_format,
            )

        folders_outfile = os.path.join(outdir_data, f"{data_set_name}_sitelist.json")
        with open(folders_outfile, "w") as file:
            json.dump(all_usefolders, file)

    printandlog("Finished 2.merge-single-cells.")


#################################
# MAIN USER INTERACTION
#################################

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Pass paths to defaults config and experiment config.")
        sys.exit()
    path_to_defaults_config = sys.argv[1]
    path_to_experiment_config = sys.argv[2]
    merge_single_cells(path_to_defaults_config, path_to_experiment_config)
