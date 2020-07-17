import os
import sys
import pathlib
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append(os.path.join("..", "scripts"))
from config_utils import process_config_file
from cell_quality_utils import CellQuality

parser = argparse.ArgumentParser()
parser.add_argument(
    "--config_file",
    help="configuration yaml file for preprocessing pipeline",
    default="site_processing_config.yaml",
)
args = parser.parse_args()
config_file = args.config_file

config = process_config_file(config_file)

# Defines the sections of the config file
core_args = config["core"]
spots_args = config["process-spots"]
summ_cell_args = config["summarize-cells"]
summ_plate_args = config["summarize-plate"]

# Defines the variables set in the config file
batch = core_args["batch"]
quality_func = core_args["categorize_cell_quality"]
ignore_files = core_args["ignore_files"]
sites_per_image_grid_side = core_args["sites_per_image_grid_side"]
batch_dir = core_args["batch_dir"]

cell_filter = spots_args["cell_filter"]

summ_cells_results_basedir = summ_cell_args["output_resultsdir"]
summ_cells_figures_basedir = summ_cell_args["output_figuresdir"]
cell_category_order = summ_cell_args["cell_category_order"]

correlation_threshold = summ_plate_args["correlation_threshold"]
painting_image_names = summ_plate_args["painting_image_names"]
barcoding_cycles = summ_plate_args["barcoding_cycles"]
barcoding_prefix = summ_plate_args["barcoding_prefix"]

cell_count_file = pathlib.Path(
    summ_cells_results_basedir, batch, "cells", "cell_count.tsv"
)
cell_count_df = pd.read_csv(cell_count_file, sep="\t")

figures_output = pathlib.Path(summ_cells_figures_basedir, batch)
os.makedirs(figures_output, exist_ok=True)
results_output = pathlib.Path(summ_cells_results_basedir, batch)
os.makedirs(results_output, exist_ok=True)

# Creates x, y coordinates for plotting per-plate views.
# Assumes image numbering starts in upper left corner and proceeds down
final_order = []
for i in range(1, sites_per_image_grid_side + 1):
    build_seq = list(
        zip(
            ([i] * (sites_per_image_grid_side + 1)),
            reversed(range(1, (sites_per_image_grid_side + 1))),
        )
    )
    final_order += build_seq

# Uses sites_list in case there are fewer analyzed sites than acquired sites
sites_list = [*range(1, (sites_per_image_grid_side * sites_per_image_grid_side) + 1)]
loc_df = (
    pd.DataFrame(final_order)
    .rename(columns={0: "x_loc", 1: "y_loc"})
    .assign(Site=sites_list)
)

# Create total_cell_count
cell_count_bysite_df = (
    cell_count_df.groupby("site_full")["cell_count"]
    .sum()
    .reset_index()
    .rename(columns={"cell_count": "total_cell_count"})
)

# Add total_cell_count to cell_count_df and add in x, y coordinates for plotting
cell_count_df = cell_count_df.merge(cell_count_bysite_df, on="site_full").merge(
    loc_df, on="Site"
)

# Plot total number of cells per well
cell_count_totalcells_df = (
    cell_count_df.groupby(["x_loc", "y_loc", "Well", "Site"])["total_cell_count"]
    .mean()
    .reset_index()
)
os.makedirs(figures_output, exist_ok=True)

by_well_gg = (
    gg.ggplot(cell_count_totalcells_df, gg.aes(x="x_loc", y="y_loc"))
    + gg.geom_point(gg.aes(fill="total_cell_count"), size=10)
    + gg.geom_text(gg.aes(label="Site"), color="lightgrey")
    + gg.facet_wrap("~Well")
    + gg.coord_fixed()
    + gg.theme_bw()
    + gg.theme(
        axis_text=gg.element_blank(),
        axis_title=gg.element_blank(),
        strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
    )
    + gg.scale_fill_cmap(name="magma")
)
output_file = pathlib.Path(figures_output, "plate_layout_cells_count_per_well.png")
by_well_gg.save(output_file, dpi=300, verbose=False)
by_well_gg

# Plot cell category ratios per well
ratio_df = pd.pivot_table(
    cell_count_df,
    values="cell_count",
    index=["site_full", "Plate", "Well", "Site", "x_loc", "y_loc"],
    columns=["Cell_Quality"],
)
ratio_df = ratio_df.assign(
    Sum=ratio_df.sum(axis=1), Pass_Filter=ratio_df[cell_filter].sum(axis=1)
)
fail_filter = [cat for cat in cell_category_order if cat not in cell_filter]
fail_filter_noempty = [
    cat for cat in cell_category_order if cat not in cell_filter if cat != "Empty"
]
ratio_df = ratio_df.assign(
    Fail_Filter=ratio_df[fail_filter].sum(axis=1),
    Fail_Filter_noempty=ratio_df[fail_filter_noempty].sum(axis=1),
)
ratio_df = ratio_df.assign(
    Pass_Fail_withempty=ratio_df["Pass_Filter"] / ratio_df["Fail_Filter"],
    Pass_Fail_0empty=ratio_df["Pass_Filter"] / ratio_df["Fail_Filter_noempty"],
    Percent_Empty=ratio_df["Empty"] / ratio_df["Sum"],
)
ratio_df = (
    ratio_df.drop(
        cell_category_order
        + ["Sum", "Pass_Filter", "Fail_Filter", "Fail_Filter_noempty"],
        1,
    )
    .stack()
    .to_frame()
    .reset_index()
    .rename(columns={0: "value"})
)
ratio_gg = (
    gg.ggplot(ratio_df, gg.aes(x="x_loc", y="y_loc"))
    + gg.geom_point(gg.aes(fill="value"), size=10)
    + gg.geom_text(gg.aes(label="Site"), color="lightgrey")
    + gg.facet_grid("Cell_Quality~Well", scales="free_y")
    + gg.coord_fixed()
    + gg.theme_bw()
    + gg.theme(
        axis_text=gg.element_blank(),
        axis_title=gg.element_blank(),
        strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
    )
    + gg.scale_fill_cmap(name="magma")
)

output_file = pathlib.Path(figures_output, "plate_layout_ratios_per_well.png")
ratio_gg.save(output_file, dpi=300, verbose=False)

# Create image dataframe
sites = [x for x in os.listdir(batch_dir) if x not in ignore_files]
image_list = []

for site in sites:
    try:
        print(f"Now processing {site}...")
        image_file = pathlib.Path(batch_dir, site, "Image.csv")
        # Aggregates image information by site into single list
        image_df = pd.read_csv(image_file)
        image_df["site"] = site
        image_df = image_df.assign(
            Plate=[x[0] for x in image_df.site.str.split("-")],
            Well=[x[1] for x in image_df.site.str.split("-")],
            Site=[x[2] for x in image_df.site.str.split("-")],
        )
        image_list.append(image_df)
    except FileNotFoundError:
        print(f"{site} data not found")
        continue

image_df = pd.concat(image_list, axis="rows").reset_index(drop=True)
print("Done concatenating image files")

# Add in x, y coordinates for plotting
image_df["Site"] = image_df["Site"].astype(int)
loc_df["Site"] = loc_df["Site"].astype(int)
image_df = image_df.merge(loc_df, how="left", on="Site")

# Plot final Cells thresholds per well
cells_finalthresh_gg = (
    gg.ggplot(image_df, gg.aes(x="x_loc", y="y_loc"))
    + gg.geom_point(gg.aes(fill="Threshold_FinalThreshold_Cells"), size=10)
    + gg.geom_text(gg.aes(label="Site"), color="lightgrey")
    + gg.facet_wrap("~Well")
    + gg.coord_fixed()
    + gg.theme_bw()
    + gg.theme(
        axis_text=gg.element_blank(),
        axis_title=gg.element_blank(),
        strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
    )
    + gg.scale_fill_cmap(name="magma")
)
output_file = pathlib.Path(
    figures_output, "plate_layout_Cells_FinalThreshold_per_well.png"
)
cells_finalthresh_gg.save(output_file, dpi=300, verbose=False)

# Plot final Nuclei thresholds per well
nuclei_finalthresh_gg = (
    gg.ggplot(image_df, gg.aes(x="x_loc", y="y_loc"))
    + gg.geom_point(gg.aes(fill="Threshold_FinalThreshold_Nuclei"), size=10)
    + gg.geom_text(gg.aes(label="Site"), color="lightgrey")
    + gg.facet_wrap("~Well")
    + gg.coord_fixed()
    + gg.theme_bw()
    + gg.theme(
        axis_text=gg.element_blank(),
        axis_title=gg.element_blank(),
        strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
    )
    + gg.scale_fill_cmap(name="magma")
)
output_file = pathlib.Path(
    figures_output, "plate_layout_Nuclei_FinalThreshold_per_well.png"
)
nuclei_finalthresh_gg.save(output_file, dpi=300, verbose=False)

# Plot percent Confluent regions per well
percent_confluent_gg = (
    gg.ggplot(image_df, gg.aes(x="x_loc", y="y_loc"))
    + gg.geom_point(gg.aes(fill="Math_PercentConfluent"), size=10)
    + gg.geom_text(gg.aes(label="Site"), color="lightgrey")
    + gg.facet_wrap("~Well")
    + gg.coord_fixed()
    + gg.theme_bw()
    + gg.theme(
        axis_text=gg.element_blank(),
        axis_title=gg.element_blank(),
        strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
    )
    + gg.scale_fill_cmap(name="magma")
)
output_file = pathlib.Path(
    figures_output, "plate_layout_PercentConfluent_per_well.png"
)
percent_confluent_gg.save(output_file, dpi=300, verbose=False)

# Create list of sites with confluent regions
confluent_df = image_df.loc[image_df["Math_PercentConfluent"] > 0]
confluent_df = (
    confluent_df[["site", "Plate", "Well", "Site", "Math_PercentConfluent"]]
    .sort_values(by=["site"])
    .reset_index(drop=True)
)
if len(confluent_df.index) > 0:
    confluent_output_file = pathlib.Path(
        results_output, "sites_with_confluent_regions.csv"
    )
    confluent_df.to_csv(confluent_output_file)

# Power Log Log Slope on Cell Painting images (proxy for focus)
PLLS_df_cols = ["Plate", "Well", "Site"]
PLLS_cols = []
for name in painting_image_names:
    PLLS_df_cols.append("ImageQuality_PowerLogLogSlope_" + name)
    PLLS_cols.append("ImageQuality_PowerLogLogSlope_" + name)
PLLS_df = image_df.loc[:, PLLS_df_cols]
PLLS_df = PLLS_df.melt(id_vars=["Plate", "Well", "Site"], var_name="Channel").replace(
    {"ImageQuality_PowerLogLogSlope_": ""}, regex=True
)

PLLS_gg = (
    gg.ggplot(PLLS_df, gg.aes(x="Site", y="value", label="Site"))
    + gg.geom_point()
    + gg.coord_fixed(ratio=0.25)
    + gg.geom_text(nudge_y=-0.025)
    + gg.facet_grid("Channel~Well", scales="free_y")
    + gg.theme_bw()
    + gg.theme(strip_background=gg.element_rect(colour="black", fill="#fdfff4"))
)

output_file = pathlib.Path(figures_output, "PLLS_per_well.png")
PLLS_gg.save(output_file, dpi=300, verbose=False)

# Outputs list of sites that are saturated in any channel
# Cell Painting images use >1% saturated, Barcoding images uses >.25% saturated
cp_sat_cols = []
bc_sat_cols = []
nts = ["A", "C", "G", "T"]
for name in painting_image_names:
    cp_sat_cols.append("ImageQuality_PercentMaximal_" + name)
for x in range(1, (barcoding_cycles + 1)):
    for nt in nts:
        bc_sat_cols.append(
            "ImageQuality_PercentMaximal_" + barcoding_prefix + "%02d" % x + "_" + nt
        )
for col in cp_sat_cols:
    cp_sat_df = image_df[image_df[col] > 1]
for col in bc_sat_cols:
    bc_sat_df = image_df[image_df[col] > 0.25]
sat_df_cols = cp_sat_cols + bc_sat_cols
sat_df_cols.append("site")

sat_df = cp_sat_df.append(bc_sat_df).drop_duplicates(subset="site")

if len(sat_df.index) > 0:
    sat_output_file = pathlib.Path(results_output, "saturated_sites.csv")
    sat_df.to_csv(sat_output_file)

# Plots saturation in Cell Painting images
cp_sat_df_cols = ["Plate", "Well", "Site"]
for name in painting_image_names:
    cp_sat_df_cols.append("ImageQuality_PercentMaximal_" + name)
    cp_sat_df_cols.append("ImageQuality_StdIntensity_" + name)
cp_sat_df = image_df.loc[:, cp_sat_df_cols]

cp_sat_df = cp_sat_df.set_index(["Plate", "Well", "Site"]).stack().reset_index()
cp_sat_df[["cat", "type", "Ch"]] = cp_sat_df["level_3"].str.split("_", n=2, expand=True)
cp_sat_df = cp_sat_df.drop(["level_3", "cat"], 1)
cp_sat_df = pd.pivot_table(
    cp_sat_df, index=["Plate", "Well", "Site", "Ch"], columns=["type"]
).reset_index()
cp_sat_df.columns = ["Plate", "Well", "Site", "Ch", "PercentMax", "StdIntensity"]

cp_saturation_gg = (
    gg.ggplot(cp_sat_df, gg.aes(x="StdIntensity", y="PercentMax", label="Site"))
    + gg.geom_point()
    + gg.coord_fixed(ratio=0.25)
    + gg.geom_text()
    + gg.facet_wrap(["Ch", "Well"], nrow=len(painting_image_names), scales="free")
    + gg.theme_bw()
    + gg.theme(strip_background=gg.element_rect(colour="black", fill="#fdfff4"))
)
output_file = pathlib.Path(figures_output, "cp_saturation.png")
cp_saturation_gg.save(output_file, dpi=300, verbose=False)

# Plots saturation in Barcoding images
bc_sat_df_cols = ["Plate", "Well", "Site"]
for x in range(1, (barcoding_cycles + 1)):
    for nt in nts:
        bc_sat_df_cols.append(
            "ImageQuality_PercentMaximal_" + barcoding_prefix + "%02d" % x + "_" + nt
        )
        bc_sat_df_cols.append(
            "ImageQuality_StdIntensity_" + barcoding_prefix + "%02d" % x + "_" + nt
        )
bc_sat_df = image_df.loc[:, bc_sat_df_cols]

bc_sat_df = bc_sat_df.set_index(["Plate", "Well", "Site"]).stack().reset_index()
bc_sat_df[["cat", "type", "Ch"]] = bc_sat_df["level_3"].str.split("_", n=2, expand=True)
bc_sat_df = bc_sat_df.drop(["level_3", "cat"], 1)
bc_sat_df = pd.pivot_table(
    bc_sat_df, index=["Plate", "Well", "Site", "Ch"], columns=["type"]
).reset_index()
bc_sat_df.columns = ["Plate", "Well", "Site", "Ch", "PercentMax", "StdIntensity"]

bc_saturation_gg = (
    gg.ggplot(bc_sat_df, gg.aes(x="StdIntensity", y="PercentMax", label="Site"))
    + gg.geom_point()
    + gg.coord_fixed(ratio=0.25)
    + gg.geom_text()
    + gg.facet_wrap(["Ch", "Well"], ncol=4, scales="free")
    + gg.theme_bw()
    + gg.theme(strip_background=gg.element_rect(colour="black", fill="#fdfff4"))
)
output_file = pathlib.Path(figures_output, "bc_saturation.png")
bc_saturation_gg.save(output_file, dpi=300, verbose=False)

# Create list of questionable channel correlations (alignments)
corr_df_cols = ["Plate", "Well", "Site", "site"]
corr_cols = []
for col in image_df.columns:
    if "Correlation_Correlation_" in col:
        corr_cols.append(col)
        corr_df_cols.append(col)
image_corr_df = image_df[corr_df_cols]
image_corr_list = []
for col in corr_cols:
    image_corr_list.append(
        image_corr_df.loc[image_corr_df[col] < correlation_threshold]
    )
image_corr_df = pd.concat(image_corr_list).drop_duplicates(subset="site").reset_index()
for col in corr_cols:
    image_corr_df.loc[(image_corr_df[col] >= correlation_threshold), col] = "pass"

if len(image_corr_df.index) > 0:
    corr_output_file = pathlib.Path(results_output, "flagged_correlations.csv")
    image_corr_df.to_csv(corr_output_file)
