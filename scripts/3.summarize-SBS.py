"""
Does not allow for selection of overwrite.
"""

import os
import sys
import logging
import traceback
import pandas as pd
import json
import math
import seaborn.objects as so

recipe_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(recipe_path, "utils"))
from io_utils import load_configs, read_csvs_with_chunksize
from cell_quality_utils import get_cell_quality_dict
from plotting_utils import make_loc_df, make_plate_layout_plots

# Configure logging
logfolder = os.path.join(recipe_path, "logs")
if not os.path.isdir(logfolder):
    os.mkdir(logfolder)
logging.basicConfig(
    filename=os.path.join(logfolder, "3.summarize-SBS.log"),
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


def make_bc_counts(barcode_list):
    df = pd.DataFrame()
    df["All_Called_Barcodes"] = barcode_list
    df = df.groupby(["All_Called_Barcodes"]).size().reset_index(name="count")
    return df


def make_summary_graph(df, id_var, categories, outpath):
    melt = df.melt(
        id_vars=id_var,
        value_vars=[x for x in df.columns if "num_quality" in x],
    )
    melt["variable"] = melt["variable"].str.replace("num_quality_", "")
    melt["variable"] = pd.Categorical(
        melt["variable"], categories=categories, ordered=True
    )
    if id_var:
        p = (
            so.Plot(melt, x="variable", y="value", color="variable")
            .facet(col=id_var, wrap=3)
            .add(so.Bar())
            .label(x="Quality Category", y="count")
        )
    else:
        p = (
            so.Plot(melt, x="variable", y="value", color="variable")
            .add(so.Bar(), so.Stack())
            .label(x="Quality Category", y="count")
        )
    p.save(outpath, dpi=300)


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


def summarize_SBS(path_to_defaults_config, path_to_experiment_config):
    printandlog("Started 3.summarize-SBS.")
    defaults_config, experiment_config = load_configs(
        path_to_defaults_config, path_to_experiment_config
    )

    # define experiment variables
    data_sets = experiment_config["data_sets"]
    cell_quality_method = experiment_config["cell_quality_method"]
    keep_cell_qualities = experiment_config["keep_cell_qualities"]

    # define defaults variables
    out_root = defaults_config["directory_structure"]["root"]
    sbs_folder = defaults_config["directory_structure"]["SBS"]
    data_folder = defaults_config["directory_structure"]["data"]
    figures_folder = defaults_config["directory_structure"]["figures"]
    cell_quality_column = defaults_config["core"]["cell_quality_column"]
    barcode_col = defaults_config["process"]["process_SBS"]["barcode_col"]
    allowed_skips = defaults_config["process"]["summarize_SBS"]["allowed_skips"]
    cell_quality_dict = get_cell_quality_dict(cell_quality_method)

    outdir_data = os.path.join(out_root, data_folder)
    if not os.path.exists(outdir_data):
        os.makedirs(outdir_data, exist_ok=True)
    outdir_figs = os.path.join(out_root, figures_folder)
    if not os.path.exists(outdir_figs):
        os.makedirs(outdir_figs, exist_ok=True)

    sitelist = []
    all_called_barcodes = []
    all_site_stats = []
    allowed_skip_counter = 0
    for data_set_name in data_sets.keys():
        per_set_called_quality = pd.DataFrame()
        printandlog(f"Starting {data_set_name}")
        folders_infile = os.path.join(outdir_data, f"{data_set_name}_sitelist.json")
        with open(folders_infile, "r") as file:
            usefolders = json.load(file)

        data_set_called_barcodes = []
        for plate_well_site_folder in usefolders:
            indir = os.path.join(out_root, sbs_folder, plate_well_site_folder)
            try:
                df = read_csvs_with_chunksize(
                    os.path.join(indir, "processed_SBS.tsv.gz"),
                    usecols=[barcode_col, "Metadata_Site"],
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
            try:
                with open(os.path.join(indir, "site_stats.json"), "r") as f:
                    site_stats = json.load(f)
            except:
                allowed_skip_counter = fail_site(
                    plate_well_site_folder,
                    "site_stats.json",
                    allowed_skip_counter,
                    allowed_skips,
                )
                continue
            try:
                indir = os.path.join(
                    out_root,
                    sbs_folder,
                    plate_well_site_folder,
                )
                quality_df = read_csvs_with_chunksize(
                    os.path.join(
                        indir, f"{plate_well_site_folder}_cell_quality_summary.csv"
                    ),
                )
            except:
                allowed_skip_counter = fail_site(
                    plate_well_site_folder,
                    f"{plate_well_site_folder}_cell_quality_summary.csv",
                    allowed_skip_counter,
                    allowed_skips,
                )
                continue

            data_set_called_barcodes += df[barcode_col].tolist()
            all_site_stats.append(site_stats)
            all_called_barcodes += data_set_called_barcodes
            per_set_called_quality = pd.concat([per_set_called_quality, quality_df])
            sitelist += [int(df["Metadata_Site"][0])]

        # Create and save barcode count summary per data set
        barcode_count_summary_df = make_bc_counts(data_set_called_barcodes)
        output_file = os.path.join(
            outdir_data, f"Total_Barcode_Calls_Counts_{data_set_name}.tsv"
        )
        barcode_count_summary_df.to_csv(output_file, sep="\t", index=False)

        # Visualize summary quality ratios per data set
        sites_per_image_grid_side = int(math.ceil(math.sqrt(max(list(set(sitelist))))))
        loc_df = make_loc_df(sites_per_image_grid_side)

        per_set_called_quality = per_set_called_quality.merge(
            loc_df, on="Metadata_Site"
        )

        ratio_df = pd.pivot_table(
            per_set_called_quality,
            values="count",
            index=[
                "Metadata_Plate",
                "Metadata_Well",
                "Metadata_Site",
                "Metadata_Dataset_Split",
                "x_loc",
                "y_loc",
            ],
            columns=[cell_quality_column],
        )
        ratio_df = ratio_df.assign(
            Sum=ratio_df.sum(axis=1),
            Pass_Filter=ratio_df[[x for x in keep_cell_qualities if x in ratio_df.columns]].sum(axis=1),
        )
        fail_filter = [
            cat
            for cat in per_set_called_quality[cell_quality_column].unique()
            if cat not in keep_cell_qualities
        ]
        fail_filter_noempty = [cat for cat in fail_filter if cat != "Empty"]
        not_empty = [
            cat
            for cat in per_set_called_quality[cell_quality_column].unique()
            if cat != "Empty"
        ]
        ratio_df = ratio_df.assign(
            Fail_Filter=ratio_df[[x for x in fail_filter if x in ratio_df.columns]].sum(axis=1),
            Fail_Filter_noempty=ratio_df[[x for x in fail_filter_noempty if x in ratio_df.columns]].sum(axis=1),
            NotEmpty=ratio_df[[x for x in not_empty if x in ratio_df.columns]].sum(axis=1),
        )
        ratio_df = ratio_df.assign(
            Pass_Fail_withempty=ratio_df["Pass_Filter"] / ratio_df["Fail_Filter"],
            Pass_Fail_0empty=ratio_df["Pass_Filter"] / ratio_df["Fail_Filter_noempty"],
            PercentEmpty=ratio_df["Empty"] / ratio_df["NotEmpty"] * 100,
        )

        try:
            outpath = os.path.join(
                outdir_figs, f"plate_layout_PassToFailRatio_{data_set_name}.png"
            )
            make_plate_layout_plots(
                ratio_df,
                "Pass_Fail_withempty",
                "Pass:Fail (with empty cells)",
                outpath,
                legend=True,
            )

            outpath = os.path.join(
                outdir_figs, f"plate_layout_PassToFailRatio_0empty_{data_set_name}.png"
            )
            make_plate_layout_plots(
                ratio_df,
                "Pass_Fail_0empty",
                "Pass:Fail (without empty cells)",
                outpath,
                legend=True,
            )

            outpath = os.path.join(
                outdir_figs, f"plate_layout_PercentEmpty_{data_set_name}.png"
            )
            make_plate_layout_plots(
                ratio_df, "PercentEmpty", "Percent Empty Cells", outpath, legend=True
            )
        except:
            printandlog(f"Failed to create Pass/Fail plots for {data_set_name}")

    # Create total barcode count summary for easy NGS comparison
    barcode_count_summary_df = make_bc_counts(all_called_barcodes)
    output_file = os.path.join(
        outdir_data, f"Total_Barcode_Calls_Counts_WholeExperiment.tsv"
    )
    barcode_count_summary_df.to_csv(output_file, sep="\t", index=False)

    # Create summary data stats
    stats_df = pd.DataFrame(all_site_stats)
    stats_df[[x for x in stats_df.columns if "Metadata" not in x]].sum().to_json(
        os.path.join(outdir_data, "Data_Stats_WholeExperiment.json"), indent=4
    )
    for data_set_name in data_sets.keys():
        df = stats_df.groupby("Metadata_Dataset_Split").sum().reset_index()
        df[[x for x in df.columns if "Metadata" not in x]].iloc[0].to_json(
            os.path.join(outdir_data, f"Data_Stats_{data_set_name}.json"),
            index=False,
            indent=4,
        )

    # Graph summary stats
    make_summary_graph(
        df,
        id_var="Metadata_Dataset_Split",
        categories=cell_quality_dict.values(),
        outpath=os.path.join(outdir_figs, "Quality_by_data_split"),
    )
    make_summary_graph(
        df,
        id_var="Metadata_Well",
        categories=cell_quality_dict.values(),
        outpath=os.path.join(outdir_figs, "Quality_by_well"),
    )
    make_summary_graph(
        df,
        id_var=[],  # All data
        categories=cell_quality_dict.values(),
        outpath=os.path.join(outdir_figs, "Quality_all_data"),
    )

    printandlog("Done with 3.summarize-SBS")


#################################
# MAIN USER INTERACTION
#################################

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Pass paths to defaults config and experiment config.")
        sys.exit()
    path_to_defaults_config = sys.argv[1]
    path_to_experiment_config = sys.argv[2]
    summarize_SBS(path_to_defaults_config, path_to_experiment_config)
