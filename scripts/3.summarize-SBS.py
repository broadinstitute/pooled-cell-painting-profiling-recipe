"""
Does not allow for selection of overwrite.

Outputs the following files:
- {figures}/plate_layout_PassToFailRatio_{data_set_name}.png
- {figures}/plate_layout_PassToFailRatio_0empty_{data_set_name}.png
- {figures}/plate_layout_PercentEmpty_{data_set_name}.png
- {figures}/Quality_Cells_by_{data_split/well/all_data}.png
- summary_data/Total_Barcode_Calls_Counts_{data_set_name}.tsv
- summary_data/Total_Barcode_Calls_Counts_WholeExperiment.tsv
- summary_data/Data_Stats_{data_set_name}.json
- summary_data/Data_Stats_WholeExperiment.json
"""

import os
import sys
import logging
import traceback
import pandas as pd
import json
import math
import seaborn.objects as so
import matplotlib.pyplot as plt

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


def make_summary_graph(df, id_var, categories, outpath, object):
    melt = df.melt(
        id_vars=id_var,
        value_vars=categories,
    )
    if 'Metadata_Quality_Name' in melt.columns:
        variable = "Metadata_Quality_Name"
    else:
        variable = "variable"
    melt[variable] = melt[variable].str.replace("num_quality_", "")
    melt[variable] = melt[variable].str.replace("_spots", "")
    categories=melt[variable].unique()
    melt[variable] = pd.Categorical(
        melt[variable], categories=categories, ordered=True
    )

    # Size the figure so each facet stays a legible fixed size no matter how
    # many facets there are, and so there's dedicated room for the legend
    # (otherwise it gets squeezed off the right edge of the figure).
    subplot_size = 3
    legend_width = 2.5
    wrap = 3
    if id_var:
        n_facets = melt[id_var].nunique()
        ncols = min(wrap, n_facets)
        nrows = math.ceil(n_facets / ncols)
        fig = plt.figure(
            figsize=(ncols * subplot_size + legend_width, nrows * subplot_size)
        )
        (
            so.Plot(melt, x=variable, y="value", color=variable)
            .facet(col=id_var, wrap=wrap)
            .add(so.Bar(alpha=1))
            .label(x="", y=f"{object} count", color=f"{object} Quality")
            .layout(engine="constrained")
            .on(fig)
            .plot()
        )
    else:
        fig = plt.figure(figsize=(subplot_size + legend_width, subplot_size))
        (
            so.Plot(melt, x=variable, y="value", color=variable)
            .add(so.Bar(alpha=1), so.Stack())
            .label(x="", y=f"{object} count", color=f"{object} Quality")
            .layout(engine="constrained")
            .on(fig)
            .plot()
        )
    for ax in fig.axes:
        ax.set_box_aspect(1)
        ax.tick_params(axis='x', labelrotation=90)
    # seaborn.objects anchors the legend using a bbox_transform that goes stale
    # during the bbox_inches="tight" save pass, which pushes it off the right
    # edge of the saved image. Re-anchoring to the figure's live transform fixes it.
    for legend in fig.legends:
        legend.set_bbox_to_anchor((0.98, 0.55), transform=fig.transFigure)
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)


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
    ratio_dict = {}
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
            Pass_Filter=ratio_df[
                [x for x in keep_cell_qualities if x in ratio_df.columns]
            ].sum(axis=1),
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
            Fail_Filter=ratio_df[[x for x in fail_filter if x in ratio_df.columns]].sum(
                axis=1
            ),
            Fail_Filter_noempty=ratio_df[
                [x for x in fail_filter_noempty if x in ratio_df.columns]
            ].sum(axis=1),
            NotEmpty=ratio_df[[x for x in not_empty if x in ratio_df.columns]].sum(
                axis=1
            ),
        )
        ratio_df = ratio_df.assign(
            Ratio_PassToFail_WithEmptyCells=ratio_df["Pass_Filter"] / ratio_df["Fail_Filter"],
            Ratio_PassToFail_WithoutEmptyCells=ratio_df["Pass_Filter"] / ratio_df["Fail_Filter_noempty"],
            PercentEmptyCells=ratio_df["Empty"] / ratio_df["Sum"] * 100,
        )

        try:
            outpath = os.path.join(
                outdir_figs, f"plate_layout_PassToFailRatio_{data_set_name}.png"
            )
            make_plate_layout_plots(
                ratio_df.reset_index(),
                "Ratio_PassToFail_WithEmptyCells",
                "Pass:Fail (with empty cells)",
                outpath,
                legend=True,
            )

            outpath = os.path.join(
                outdir_figs, f"plate_layout_PassToFailRatio_0empty_{data_set_name}.png"
            )
            make_plate_layout_plots(
                ratio_df.reset_index(),
                "Ratio_PassToFail_WithoutEmptyCells",
                "Pass:Fail (without empty cells)",
                outpath,
                legend=True,
            )

            outpath = os.path.join(
                outdir_figs, f"plate_layout_PercentEmpty_{data_set_name}.png"
            )
            make_plate_layout_plots(
                ratio_df.reset_index(),
                "PercentEmptyCells",
                "Percent Empty Cells",
                outpath,
                legend=True,
            )
        except:
            printandlog(f"Failed to create Pass/Fail plots for {data_set_name}")
        
        ratio_dict[data_set_name] = ratio_df

    ratio_df = pd.concat(ratio_dict)

    # Create total barcode count summary for easy NGS comparison
    barcode_count_summary_df = make_bc_counts(all_called_barcodes)
    output_file = os.path.join(
        outdir_data, f"Total_Barcode_Calls_Counts_WholeExperiment.tsv"
    )
    barcode_count_summary_df.to_csv(output_file, sep="\t", index=False)

    stats_df = pd.DataFrame(all_site_stats)

    # Plot summary barcode stats
    make_summary_graph(
        stats_df[[x for x in stats_df.columns if 'num_quality' in x and 'spots' in x]+["Metadata_Dataset_Split"]],
        id_var="Metadata_Dataset_Split",
        categories=[x for x in stats_df.columns if 'num_quality' in x and 'spots' in x],
        outpath=os.path.join(outdir_figs, "Quality_Spots_by_data_split"),
        object="Foci"
    )
    make_summary_graph(
        stats_df[[x for x in stats_df.columns if 'num_quality' in x and 'spots' in x]+["Metadata_Well"]],
        id_var="Metadata_Well",
        categories=[x for x in stats_df.columns if 'num_quality' in x and 'spots' in x],
        outpath=os.path.join(outdir_figs, "Quality_Spots_by_well"),
        object="Foci"
    )
    make_summary_graph(
        stats_df[[x for x in stats_df.columns if 'num_quality' in x and 'spots' in x]],
        id_var=[],  # All data
        categories=[x for x in stats_df.columns if 'num_quality' in x and 'spots' in x],
        outpath=os.path.join(outdir_figs, "Quality_Spots_all_data"),
        object="Foci"
    )

    # Create summary data stats
    ratio_df = ratio_df.reset_index()
    for data_set_name in data_sets.keys():
        df_sum = stats_df.loc[stats_df['Metadata_Dataset_Split']==data_set_name][[x for x in stats_df.columns if 'percent' not in x]].groupby("Metadata_Dataset_Split").sum().reset_index()
        df_pct = stats_df.loc[stats_df['Metadata_Dataset_Split']==data_set_name][[x for x in stats_df.columns if 'percent' in x]+["Metadata_Dataset_Split"]].groupby("Metadata_Dataset_Split").mean().reset_index()
        df = pd.concat([df_sum[[x for x in df_sum.columns if "Metadata" not in x]], df_pct[[x for x in df_pct.columns if "Metadata" not in x]]], axis=1).iloc[0]
        df['num_empty_cells'] = ratio_df.loc[ratio_df['Metadata_Dataset_Split']==data_set_name]['Empty'].sum()
        df['num_not_empty_cells'] = ratio_df.loc[ratio_df['Metadata_Dataset_Split']==data_set_name]['NotEmpty'].sum()
        df['num_total_cells'] = ratio_df.loc[ratio_df['Metadata_Dataset_Split']==data_set_name]['Sum'].sum()
        df['percent_empty_cells'] = df['num_empty_cells'] / df['num_total_cells']
        # remove data that shouldn't be summed or averaged
        df.pop("num_unique_genes")
        df.pop("num_unique_guides")
        df.to_json(
            os.path.join(outdir_data, f"Data_Stats_{data_set_name}.json"),
            index=False,
            indent=4,
        )
    df_sum = stats_df[[x for x in stats_df.columns if 'percent' not in x if 'Metadata' not in x]].sum()
    df_pct = stats_df[[x for x in stats_df.columns if 'percent' in x if 'Metadata' not in x]].mean()
    stats_df = pd.concat([df_sum, df_pct])
    stats_df['num_empty_cells'] = ratio_df['Empty'].sum()
    stats_df['num_not_empty_cells'] = ratio_df['NotEmpty'].sum()
    stats_df['num_total_cells'] = ratio_df['Sum'].sum()
    stats_df['percent_empty_cells'] = stats_df['num_empty_cells'] / stats_df['num_total_cells']
    stats_df["num_spots_per_cell"] = stats_df["num_spots_in_cells"]/stats_df["num_not_empty_cells"]
    stats_df["num_unique_guides_called"] = barcode_count_summary_df.loc[barcode_count_summary_df['All_Called_Barcodes'] != "Unmatched","count"].sum()
    # remove data that shouldn't be summed or averaged
    stats_df.pop("num_unique_genes")
    stats_df.pop("num_unique_guides")
    stats_df.to_json(
        os.path.join(outdir_data, "Data_Stats_WholeExperiment.json"), indent=4
    )

    # Graph cell quality summary stats
    make_summary_graph(
        ratio_df[[x for x in cell_quality_dict.values()]+["Metadata_Dataset_Split"]],
        id_var="Metadata_Dataset_Split",
        categories=cell_quality_dict.values(),
        outpath=os.path.join(outdir_figs, "Quality_Cells_by_data_split"),
        object="Cell"
    )
    make_summary_graph(
        ratio_df[[x for x in cell_quality_dict.values()]+["Metadata_Well"]],
        id_var="Metadata_Well",
        categories=cell_quality_dict.values(),
        outpath=os.path.join(outdir_figs, "Quality_Cells_by_well"),
        object="Cell"
    )
    make_summary_graph(
        ratio_df[[x for x in cell_quality_dict.values()]],
        id_var=[],  # All data
        categories=cell_quality_dict.values(),
        outpath=os.path.join(outdir_figs, "Quality_Cells_all_data"),
        object="Cell"
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
