"""
Allows for resuming and either overwriting or not.
"""

import os
import sys
import logging
import traceback
import pandas as pd
import json

recipe_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(recipe_path, "utils"))
import barcode_calling_utils
from io_utils import load_configs, parse_data_set, read_csvs_with_chunksize
from cell_quality_utils import CellQuality


# Configure logging
logfolder = os.path.join(recipe_path, "logs")
if not os.path.isdir(logfolder):
    os.mkdir(logfolder)
logging.basicConfig(
    filename=os.path.join(logfolder, "1.process-SBS.log"),
    level=logging.INFO,
)


def printandlog(msg, type="info"):
    print(msg)
    if type == "warning":
        logging.warning(msg)
    else:
        logging.info(msg)


def handle_excepthook(exc_type, exc_value, exc_traceback):
    logging.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))
    traceback_details = "\n".join(traceback.extract_tb(exc_traceback).format())
    print(f"Uncaught Exception: {traceback_details}")


sys.excepthook = handle_excepthook


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
            f"Allowed skip limit of {allowed_skips} reached. Stopping 1.process-SBS.",
            type="warning",
        )
        sys.exit(1)
    else:
        return allowed_skip_counter


def process_SBS(path_to_defaults_config, path_to_experiment_config):
    printandlog("Started 1.process-SBS.")
    defaults_config, experiment_config = load_configs(
        path_to_defaults_config, path_to_experiment_config
    )

    # define experiment variables
    file_location = experiment_config["file_location"]
    drop_barcodes = experiment_config["drop_barcodes"]
    control_genes = experiment_config["control_genes"]
    data_sets = experiment_config["data_sets"]
    overwrite_files = experiment_config["overwrite_files"]
    compartments = experiment_config["compartments"]
    cell_quality_method = experiment_config["cell_quality_method"]
    match_to_library = experiment_config["match_to_library"]
    library_structure = experiment_config["library_structure"]
    library_location = experiment_config["library_location"]

    # define defaults variables
    out_root = defaults_config["directory_structure"]["root"]
    sbs_folder = defaults_config["directory_structure"]["SBS"]
    data_folder = defaults_config["directory_structure"]["data"]
    cell_quality_column = defaults_config["core"]["cell_quality_column"]
    cell_quality_index = defaults_config["core"]["cell_quality_index"]
    id_cols = defaults_config["core"]["cell_id_cols"]
    location_cols = defaults_config["process"]["process_SBS"]["location_cols"]
    allowed_skips = defaults_config["process"]["process_SBS"]["allowed_skips"]
    foci_cols = defaults_config["process"]["process_SBS"]["foci_cols"]
    SBS_score_col = defaults_config["process"]["process_SBS"]["SBS_score_col"]
    barcode_col = defaults_config["process"]["process_SBS"]["barcode_col"]
    gene_col = defaults_config["process"]["process_SBS"]["gene_col"]
    call_col = defaults_config["process"]["process_SBS"]["call_col"]
    spot_quality_method = defaults_config["process"]["process_SBS"]["spot_quality_method"]

    cell_quality = CellQuality(
        cell_quality_method,
        category_col_index=cell_quality_index,
        category_class_name=cell_quality_column,
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
        no_SBS_foci_sites = []
        no_assigned_cells_sites = []
        for plate in plate_list:
            for plate_well_tuple in [x for x in list_plate_well_tuples if plate in x]:
                well = plate_well_tuple[1]
                usefolders = [x for x in folderlist if f"{plate}-{well}" in x]
                usefolders = [x for x in usefolders if x not in inferred_empty_sites]

                printandlog(f"Now processing {data_set_name}: {batch}, {plate}, {well}")
                printandlog(f"{len(usefolders)} site folders to process")

                for plate_well_site_folder in usefolders:
                    output_dir = os.path.join(
                        out_root, sbs_folder, plate_well_site_folder
                    )
                    site = plate_well_site_folder.rsplit("-", 1)[1]

                    if not overwrite_files:
                        if os.path.exists(output_dir):
                            printandlog(
                                f"Skipping {plate_well_site_folder}. Output folder {output_dir} exists and overwrite_files is False."
                            )
                            continue

                    # SBS DATA HANDLING
                    try:
                        if not match_to_library:
                            barcode_file = os.path.join(
                                file_location,
                                batch,
                                plate_well_site_folder,
                                "BarcodeFoci.csv",
                            )
                            barcodefoci_df = read_csvs_with_chunksize(barcode_file)
                            # most columns from BarcodeFoci.csv are dropped
                            barcodefoci_df = barcodefoci_df[
                                id_cols
                                + location_cols
                                + [x for x in barcodefoci_df.columns if "Parent" in x]
                                + [x for x in barcodefoci_df.columns if "Child" in x]
                            ]
                        foci_file = os.path.join(
                            file_location, batch, plate_well_site_folder, "Foci.csv"
                        )
                        foci_df = read_csvs_with_chunksize(foci_file)

                        if len(foci_df) == 0:
                            printandlog(f"{plate_well_site_folder} does not have any foci")
                            no_SBS_foci_sites.append(plate_well_site_folder)
                            continue

                        if match_to_library:
                            matchdict = barcode_calling_utils.match_barcode_to_library(library_location, library_structure, call_col, foci_df)
                            for col_to_match in library_structure.keys():
                                if len(library_structure) > 1:
                                    foci_df[f"{SBS_score_col}_{col_to_match}"] = matchdict[col_to_match][0]
                                    foci_df[f"{barcode_col}_{col_to_match}"] = matchdict[col_to_match][1]
                                    foci_df[f"{gene_col}_{col_to_match}"] = matchdict[col_to_match][2]
                                    foci_df[f"{foci_cols[1]}_{col_to_match}"] = matchdict[col_to_match][3]
                                else:
                                    foci_df[SBS_score_col] = matchdict[col_to_match][0]
                                    foci_df[barcode_col] = matchdict[col_to_match][1]
                                    foci_df[gene_col] = matchdict[col_to_match][2]
                                    foci_df[foci_cols[1]] = matchdict[col_to_match][3]
                        # most columns from Foci.csv are dropped
                        foci_df = foci_df[
                            list(
                                set(
                                    id_cols
                                    + location_cols
                                    + [x for x in foci_df.columns if any(y in x for y in foci_cols)]
                                    + [x for x in foci_df.columns if any(y in x for y in [SBS_score_col, barcode_col, gene_col])]
                                )
                            )
                            + [x for x in foci_df.columns if "Parent" in x]
                            + [x for x in foci_df.columns if "Child" in x]
                        ]

                        # get length of barcode calls and confirm they are all the same length
                        if foci_df[call_col].astype(str).str.len().nunique() == 1:
                            SBScycles = foci_df[call_col].str.len()[0]
                        else:
                            printandlog(
                                f"Failed to parse number SBS cycles in {plate_well_site_folder}",
                                type="warning",
                            )
                            allowed_skip_counter = fail_site(
                                plate_well_site_folder,
                                barcode_file,
                                allowed_skip_counter,
                                allowed_skips,
                            )
                            continue
                        
                        # Add foci quality categories. Used for barcode calling if multiple barcodes. Detects recombination.
                        if len(library_structure) > 1:
                            foci_df = barcode_calling_utils.categorize_spots(foci_df, [x for x in foci_df.columns if SBS_score_col in x], SBScycles, spot_quality_method)

                    except:
                        allowed_skip_counter = fail_site(
                            plate_well_site_folder,
                            barcode_file,
                            allowed_skip_counter,
                            allowed_skips,
                        )
                        continue

                    if not match_to_library:
                        try:
                            # Confirm that image number and object number are aligned
                            pd.testing.assert_frame_equal(
                                barcodefoci_df.loc[:, id_cols],
                                foci_df.loc[:, id_cols],
                                check_names=True,
                            )
                            # Confirm that X and Y locations are aligned
                            pd.testing.assert_frame_equal(
                                barcodefoci_df.loc[:, location_cols],
                                foci_df.loc[:, location_cols],
                                check_names=True,
                            )
                        except AssertionError:
                            allowed_skip_counter = fail_site(
                                plate_well_site_folder,
                                barcode_file,
                                allowed_skip_counter,
                                allowed_skips,
                            )
                            continue

                        image_number = foci_df.ImageNumber.unique()[0]

                        # Merge SBS files
                        complete_foci_df = barcodefoci_df.merge(
                            foci_df,
                            left_on=id_cols + location_cols,
                            right_on=id_cols + location_cols,
                            how="inner",
                        )
                    else:
                        complete_foci_df = foci_df.copy()
                        image_number = complete_foci_df.ImageNumber.unique()[0]

                    if len(library_structure) == 1:
                        # Drop foci from droplist
                        complete_foci_df = complete_foci_df.loc[
                            ~complete_foci_df[barcode_col].isin(drop_barcodes)
                        ]
                    else:
                        # TODO - support dropping specific barcodes with multi-matches
                        printandlog("Dropping specific barcodes not currently supported for multiple library segments. Skipping barcode dropping.")

                    # Count foci outside of parent compartment (e.g. Cells)
                    try:
                        unassigned_spot_df = complete_foci_df.loc[
                            (
                                complete_foci_df.loc[:, f"Parent_{compartments[0]}"]
                                == 0
                            ).squeeze(),
                            :,
                        ]
                        num_unassigned_spots = unassigned_spot_df.shape[0]
                    except:
                        num_unassigned_spots = 0

                    # Count foci in parent compartment (e.g. Cells)
                    try:
                        assigned_spot_df = complete_foci_df.loc[
                            (
                                complete_foci_df.loc[:, f"Parent_{compartments[0]}"]
                                != 0
                            ).squeeze(),
                            :,
                        ]
                    except:
                        printandlog(
                            f"{plate_well_site_folder} has no foci in cells. Skipping."
                        )
                        continue

                    num_assigned_cells = len(set(
                        assigned_spot_df.loc[:, f"Parent_{compartments[0]}"]
                        )
                    )

                    if num_assigned_cells == 0:
                        printandlog(
                            f"{plate_well_site_folder} has no assigned cells. Skipping."
                        )
                        no_assigned_cells_sites.append(plate_well_site_folder)
                        continue

                    num_assigned_spots = assigned_spot_df.shape[0]

                    # Assign Cell Quality scores based on gene and barcode assignments
                    crispr_barcode_gene_df = cell_quality.assign_cell_quality(
                        assigned_spot_df,
                        parent_col=f"Parent_{compartments[0]}",
                        score_col=SBS_score_col,
                        gene_col = gene_col,
                        barcode_col = barcode_col,
                        SBScycles=SBScycles,
                        match_to_library=match_to_library,
                        library_structure=library_structure,
                    ).assign(
                        Metadata_ImageNumber=image_number,
                        Metadata_Batch=batch,
                        Metadata_Plate=plate,
                        Metadata_Well=well,
                        Metadata_Site=site,
                        Metadata_Identifier=f"{batch}-{plate}-{well}-{site}",
                        Metadata_Dataset_Split=data_set_name,
                    )

                    # SAVE OUT PROCESSED SBS DATA
                    if not os.path.exists(output_dir):
                        os.makedirs(output_dir, exist_ok=True)
                    # Used in subsequent steps to make profiles
                    out_file = os.path.join(output_dir, "processed_SBS.tsv.gz")
                    crispr_barcode_gene_df.to_csv(
                        out_file, sep="\t", index=False, compression="gzip"
                    )

                    # Create additional data summaries for subsequent qc
                    num_unique_guides = len(
                        set(crispr_barcode_gene_df.loc[:, barcode_col].to_list())
                    )
                    num_unique_genes = len(
                        set(crispr_barcode_gene_df.loc[:, gene_col].to_list())
                    )
                    gene_category_count_df = (
                        cell_quality.summarize_perturbation_quality_counts(
                            quality_df=crispr_barcode_gene_df,
                            parent_col=f"Parent_{compartments[0]}",
                            group_cols=[gene_col],
                        )
                    )

                    guide_category_count_df = (
                        cell_quality.summarize_perturbation_quality_counts(
                            quality_df=crispr_barcode_gene_df,
                            parent_col=f"Parent_{compartments[0]}",
                            group_cols=[gene_col, barcode_col],
                            guide=True,
                        )
                    )

                    count_merge_cols = list(
                        set(gene_category_count_df.columns).intersection(
                            guide_category_count_df.columns
                        )
                    )

                    cell_category_counts_df = guide_category_count_df.merge(
                        gene_category_count_df, on=count_merge_cols, how="left"
                    ).assign(
                        Metadata_ImageNumber=image_number,
                        Metadata_Batch=batch,
                        Metadata_Plate=plate,
                        Metadata_Well=well,
                        Metadata_Site=site,
                        Metadata_Identifier=f"{batch}-{plate}-{well}-{site}",
                        Metadata_Dataset_Split=data_set_name,
                    )

                    out_file = os.path.join(
                        output_dir, "cell_perturbation_category_summary_counts.tsv"
                    )
                    cell_category_counts_df.to_csv(out_file, sep="\t", index=False)

                    passed_gene_df = (
                        gene_category_count_df.groupby(gene_col)["Cell_Count_Per_Gene"]
                        .sum()
                        .reset_index()
                        .sort_values(by="Cell_Count_Per_Gene", ascending=False)
                        .reset_index(drop=True)
                    )

                    passed_gene_df.loc[:, gene_col] = pd.Categorical(
                        passed_gene_df.loc[:, gene_col].to_list(),
                        categories=passed_gene_df.loc[:, gene_col].to_list(),
                    )

                    # Number of non-targetting controls
                    nt_gene_df = passed_gene_df.query(f"{gene_col} in @control_genes")
                    num_nt = nt_gene_df.Cell_Count_Per_Gene.sum()

                    # Output descriptive stats/site
                    descriptive_results = {
                        "Metadata_ImageNumber": int(image_number),
                        "Metadata_Batch": batch,
                        "Metadata_Plate": plate,
                        "Metadata_Well": well,
                        "Metadata_Site": int(site),
                        "Metadata_Identifier": f"{batch}-{plate}-{well}-{site}",
                        "Metadata_Dataset_Split": data_set_name,
                        "num_unassigned_spots": int(num_unassigned_spots),
                        "num_assigned_spots": int(num_assigned_spots),
                        "num_unique_genes": int(num_unique_genes),
                        "num_unique_guides": int(num_unique_guides),
                        "num_assigned_cells": int(num_assigned_cells),
                        "num_nontarget_controls_kept_cells": int(num_nt),
                    }

                    cell_quality_summary_df = (
                        cell_quality.summarize_cell_quality_counts(
                            quality_df=crispr_barcode_gene_df,
                            parent_col=f"Parent_{compartments[0]}",
                        )
                    )

                    for quality in cell_quality_summary_df[cell_quality_column]:
                        descriptive_results[f"num_quality_{quality}"] = int(
                            cell_quality_summary_df.loc[
                                cell_quality_summary_df[cell_quality_column] == quality,
                                "Cell_Count",
                            ].squeeze()
                        )

                    output_file = os.path.join(output_dir, "site_stats.json")
                    with open(output_file, "w") as f:
                        json.dump(descriptive_results, f, indent=4)

        if no_assigned_cells_sites:
            no_assigned_cells_sites = pd.DataFrame(
                no_assigned_cells_sites, columns=["No_Assigned_Cells_Site_Folder"]
            )
            no_assigned_cells_sites.to_csv(
                os.path.join(
                    out_root,
                    data_folder,
                    f"No_Assigned_Cells_Sites_{data_set_name}.csv",
                ),
                index=False,
            )
        if no_SBS_foci_sites:
            no_SBS_foci_sites = pd.DataFrame(
                no_SBS_foci_sites, columns=["No_SBS_Foci_Site_Folder"]
            )
            no_SBS_foci_sites.to_csv(
                os.path.join(
                    out_root,
                    data_folder,
                    f"No_SBS_Foci_Sites_{data_set_name}.csv",
                ),
                index=False,
            )
    printandlog(f"Finished 1.process-SBS.")


#################################
# MAIN USER INTERACTION
#################################

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Pass paths to defaults config and experiment config.")
        sys.exit()
    path_to_defaults_config = sys.argv[1]
    path_to_experiment_config = sys.argv[2]
    process_SBS(path_to_defaults_config, path_to_experiment_config)
