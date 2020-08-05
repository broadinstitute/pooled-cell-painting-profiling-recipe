import os
import sys
import pathlib
import argparse
import pandas as pd
import plotnine as gg

sys.path.append("config")
from config_utils import process_config_file

recipe_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(recipe_path, "scripts"))
from cell_quality_utils import CellQuality
from arg_utils import parse_command_args
from io_utils import check_if_write

args = parse_command_args(config_file="site_processing_config.yaml")
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
image_cols = spots_args["image_cols"]
input_image_file = spots_args["image_file"]

summ_cells_results_basedir = summ_cell_args["output_summary_resultsdir"]
summ_cells_figures_basedir = summ_cell_args["output_summary_figuresdir"]
cell_category_order = summ_cell_args["cell_category_order"]
cell_count_file = summ_cell_args["cell_count_file"]

correlation_threshold = summ_plate_args["correlation_threshold"]
painting_image_names = summ_plate_args["painting_image_names"]
barcoding_cycles = summ_plate_args["barcoding_cycles"]
barcoding_prefix = summ_plate_args["barcoding_prefix"]
force = summ_plate_args["force_overwrite"]

# Forced overwrite can be achieved in one of two ways.
# The command line overrides the config file, check here if it is provided
if not force:
    force = args.force

cell_count_df = pd.read_csv(cell_count_file, sep="\t")

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
    .assign(site=sites_list)
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
    loc_df, on="site"
)

# Plot total number of cells per well
cell_count_totalcells_df = (
    cell_count_df.groupby(["x_loc", "y_loc", "well", "site"])["total_cell_count"]
    .mean()
    .reset_index()
)
os.makedirs(figures_output, exist_ok=True)

plate = cell_count_df["plate"].unique()[0]

by_well_gg = (
    gg.ggplot(cell_count_totalcells_df, gg.aes(x="x_loc", y="y_loc"))
    + gg.geom_point(gg.aes(fill="total_cell_count"), size=10)
    + gg.geom_text(gg.aes(label="site"), color="lightgrey")
    + gg.facet_wrap("~well")
    + gg.coord_fixed()
    + gg.theme_bw()
    + gg.ggtitle(f"Total Cells/Well\n{plate}")
    + gg.theme(
        axis_text=gg.element_blank(),
        axis_title=gg.element_blank(),
        strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
    )
    + gg.labs(fill="Cells")
    + gg.scale_fill_cmap(name="magma")
)

os.makedirs(figures_output, exist_ok=True)
output_file = pathlib.Path(figures_output, "plate_layout_cells_count_per_well.png")
if check_if_write(output_file, force, throw_warning=True):
    by_well_gg.save(output_file, dpi=300, verbose=False)

# Plot cell category ratios per well
ratio_df = pd.pivot_table(
    cell_count_df,
    values="cell_count",
    index=["site_full", "plate", "well", "site", "x_loc", "y_loc"],
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
    .rename(columns={0: "Ratio"})
)

quality_recode = {
    "Pass_Fail_withempty": "Pass/Fail (with empty)",
    "Pass_Fail_0empty": "Pass/Fail (without empty)",
    "Percent_Empty": "Percent Empty",
}
ratio_df = ratio_df.assign(
    cell_quality_recode=ratio_df.Cell_Quality.replace(quality_recode)
)

ratio_gg = (
    gg.ggplot(ratio_df, gg.aes(x="x_loc", y="y_loc"))
    + gg.geom_point(gg.aes(fill="Ratio"), size=5)
    + gg.geom_text(gg.aes(label="site"), size=4, color="lightgrey")
    + gg.facet_grid("cell_quality_recode~well", scales="free_y")
    + gg.coord_fixed()
    + gg.theme_bw()
    + gg.ggtitle(f"Quality Ratio\n{plate}")
    + gg.coord_fixed()
    + gg.theme(
        axis_text=gg.element_blank(),
        axis_title=gg.element_blank(),
        strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
    )
    + gg.scale_fill_cmap(name="magma")
)

output_file = pathlib.Path(figures_output, "plate_layout_ratios_per_well.png")
if check_if_write(output_file, force, throw_warning=True):
    ratio_gg.save(
        output_file,
        dpi=300,
        width=(len(ratio_df["well"].unique()) + 2),
        height=6,
        verbose=False,
    )

# Load image file
image_df = pd.read_csv(input_image_file, sep="\t")
image_meta_col_list = list(image_cols.values())

# Add in x, y coordinates for plotting
image_df["site"] = image_df[image_cols["site"]].astype(int)
loc_df["site"] = loc_df["site"].astype(int)
image_df = image_df.merge(loc_df, how="left", on="site")

# Plot final compartment thresholds per well
threshold_col_prefix = "Threshold_FinalThreshold_"
for threshhold_compartment in ["Cells", "Nuclei"]:
    threshold_full_col = f"{threshold_col_prefix}{threshhold_compartment}"

    if threshold_full_col not in image_df.columns:
        continue

    compartment_finalthresh_gg = (
        gg.ggplot(image_df, gg.aes(x="x_loc", y="y_loc"))
        + gg.geom_point(gg.aes(fill=f"{threshold_full_col}"), size=10)
        + gg.geom_text(gg.aes(label="site"), color="lightgrey")
        + gg.facet_wrap(f"~{image_cols['well']}")
        + gg.coord_fixed()
        + gg.theme_bw()
        + gg.ggtitle(f"{threshhold_compartment} Thresholds \n {plate}")
        + gg.theme(
            axis_text=gg.element_blank(),
            axis_title=gg.element_blank(),
            strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
        )
        + gg.scale_fill_cmap(name="magma")
        + gg.labs(fill="Threshold")
    )
    output_file = pathlib.Path(
        figures_output,
        f"plate_layout_{threshhold_compartment}_FinalThreshold_per_well.png",
    )
    if check_if_write(output_file, force, throw_warning=True):
        compartment_finalthresh_gg.save(output_file, dpi=300, verbose=False)

# Create list of sites with confluent regions
confluent_col = "Math_PercentConfluent"
if confluent_col in image_df.columns:
    # Plot percent confluent regions per well
    percent_confluent_gg = (
        gg.ggplot(image_df, gg.aes(x="x_loc", y="y_loc"))
        + gg.geom_point(gg.aes(fill=confluent_col), size=10)
        + gg.geom_text(gg.aes(label="site"), color="lightgrey")
        + gg.facet_wrap(f"~{image_cols['well']}")
        + gg.coord_fixed()
        + gg.theme_bw()
        + gg.ggtitle(f"Percent Confluent \n {plate}")
        + gg.theme(
            axis_text=gg.element_blank(),
            axis_title=gg.element_blank(),
            strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
        )
        + gg.scale_fill_cmap(name="magma")
        + gg.labs(fill="Percent")
    )
    output_file = pathlib.Path(
        figures_output, "plate_layout_PercentConfluent_per_well.png"
    )
    if check_if_write(output_file, force, throw_warning=True):
        percent_confluent_gg.save(output_file, dpi=300, verbose=False)

    confluent_cols = ["site", confluent_col] + image_meta_col_list
    confluent_df = image_df.loc[image_df[confluent_col] > 0]
    confluent_df = (
        confluent_df[confluent_cols].sort_values(by=["site"]).reset_index(drop=True)
    )

    if len(confluent_df.index) > 0:
        confluent_output_file = pathlib.Path(
            results_output, "sites_with_confluent_regions.csv"
        )
        if check_if_write(confluent_output_file, force, throw_warning=True):
            confluent_df.to_csv(confluent_output_file)

# Power Log Log Slope on Cell Painting images (proxy for focus)
# Any point too high or too low may have focus issues
pll_col_prefix = "ImageQuality_PowerLogLogSlope_"
PLLS_df_cols = image_meta_col_list.copy()
PLLS_cols = []
for name in painting_image_names:
    PLLS_df_cols.append(f"{pll_col_prefix}" + name)
    PLLS_cols.append(f"{pll_col_prefix}" + name)

if all(x in image_df.columns.tolist() for x in PLLS_df_cols):
    PLLS_df = image_df.loc[:, PLLS_df_cols]
    PLLS_df = PLLS_df.melt(id_vars=image_meta_col_list, var_name="channel").replace(
        {pll_col_prefix: ""}, regex=True
    )

    PLLS_gg = (
        gg.ggplot(
            PLLS_df, gg.aes(x=image_cols["site"], y="value", label=image_cols["site"])
        )
        + gg.coord_fixed(ratio=0.25)
        + gg.geom_text(size=6)
        + gg.facet_grid(f"channel~{image_cols['well']}", scales="free_y")
        + gg.ggtitle(f"Image focus \n {plate}")
        + gg.theme_bw()
        + gg.ylab("Power log log slope")
        + gg.theme(
            strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
            axis_text_y=gg.element_text(size=4),
        )
    )

    output_file = pathlib.Path(figures_output, "PLLS_per_well.png")
    if check_if_write(output_file, force, throw_warning=True):
        PLLS_gg.save(
            output_file,
            width=(len(PLLS_df[image_cols["well"]].unique()) + 2),
            height=8,
            dpi=300,
            verbose=False,
        )

# Outputs list of sites that are saturated in any channel
# Cell Painting images use >1% saturated, Barcoding images uses >.25% saturated
saturated_col_prefix = "ImageQuality_PercentMaximal_"
cp_sat_cols = []
bc_sat_cols = []
nts = ["A", "C", "G", "T"]

for name in painting_image_names:
    cp_sat_cols.append(f"{saturated_col_prefix}{name}")

for x in range(1, (barcoding_cycles + 1)):
    for nt in nts:
        bc_sat_cols.append(f"{saturated_col_prefix}{barcoding_prefix}{x:02d}_{nt}")

if all(x in image_df.columns.tolist() for x in cp_sat_cols):
    for col in cp_sat_cols:
        cp_sat_df = image_df[image_df[col] > 1]

    for col in cp_sat_cols:
        cp_sat_df = image_df[image_df[col] > 1]
    for col in bc_sat_cols:
        bc_sat_df = image_df[image_df[col] > 0.25]
    sat_df_cols = cp_sat_cols + bc_sat_cols
    sat_df_cols.append("site")

    sat_df = cp_sat_df.append(bc_sat_df).drop_duplicates(subset="site")

    if len(sat_df.index) > 0:
        sat_output_file = pathlib.Path(results_output, "saturated_sites.csv")
        if check_if_write(output_file, force, throw_warning=True):
            sat_df.to_csv(sat_output_file)

# Plots saturation in Cell Painting images
# x = std dev of intensity (to find images that have unusually bright spots)
# y = % image that is saturated (to find images that are unusually bright)
# Look at points off cluster where y > 1
intensity_col_prefix = "ImageQuality_StdIntensity_"
cp_sat_df_cols = image_meta_col_list.copy()
for name in painting_image_names:
    cp_sat_df_cols.append(f"{saturated_col_prefix}{name}")
    cp_sat_df_cols.append(f"{intensity_col_prefix}{name}")

if all(x in image_df.columns.tolist() for x in cp_sat_df_cols):
    cp_sat_df = image_df.loc[:, cp_sat_df_cols]

    cp_sat_df = cp_sat_df.set_index(image_meta_col_list).stack().reset_index()
    cp_sat_df[["cat", "type", "Ch"]] = cp_sat_df["level_3"].str.split(
        "_", n=2, expand=True
    )
    cp_sat_df = cp_sat_df.drop(["level_3", "cat"], 1)
    cp_sat_df = pd.pivot_table(
        cp_sat_df, index=image_meta_col_list + ["Ch"], columns=["type"]
    ).reset_index()
    cp_sat_df.columns = image_meta_col_list + ["Ch", "PercentMax", "StdIntensity"]

    cp_saturation_ymax = max(cp_sat_df.PercentMax)
    if cp_saturation_ymax < 1:
        cp_saturation_ymax = 1

    cp_saturation_gg = (
        gg.ggplot(
            cp_sat_df,
            gg.aes(x="StdIntensity", y="PercentMax", label=image_cols["site"]),
        )
        + gg.coord_fixed(ratio=0.25)
        + gg.geom_text(size=6)
        + gg.ylim([0, cp_saturation_ymax])
        + gg.facet_wrap(
            ["Ch", image_cols["well"]], nrow=len(painting_image_names), scales="free"
        )
        + gg.theme_bw()
        + gg.ggtitle(f"Cell Painting Image Saturation \n {plate}")
        + gg.theme(
            strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
            strip_text=gg.element_text(size=7),
            axis_text=gg.element_text(size=6),
            subplots_adjust={"wspace": 0.2},
        )
    )
    output_file = pathlib.Path(figures_output, "cp_saturation.png")
    if check_if_write(output_file, force, throw_warning=True):
        cp_saturation_gg.save(
            output_file,
            dpi=300,
            width=(len(cp_sat_df[image_cols["well"]].unique()) + 2),
            height=(len(cp_sat_df["Ch"].unique())),
            verbose=False,
        )
# Plots saturation in Barcoding images
# x = std dev of intensity (to find images that have unusually bright spots)
# y = % image that is saturated (to find images that are unusually bright)
# Look at points off cluster where y > .2
bc_sat_df_cols = image_meta_col_list.copy()
for x in range(1, (barcoding_cycles + 1)):
    for nt in nts:
        bc_sat_df_cols.append(f"{saturated_col_prefix}{barcoding_prefix}{x:02d}_{nt}")
        bc_sat_df_cols.append(f"{intensity_col_prefix}{barcoding_prefix}{x:02d}_{nt}")

if all(x in image_df.columns.tolist() for x in bc_sat_df_cols):
    bc_sat_df = image_df.loc[:, bc_sat_df_cols]

    bc_sat_df = bc_sat_df.set_index(image_meta_col_list).stack().reset_index()
    bc_sat_df[["cat", "type", "Ch"]] = bc_sat_df["level_3"].str.split(
        "_", n=2, expand=True
    )
    bc_sat_df = bc_sat_df.drop(["level_3", "cat"], 1)
    bc_sat_df = pd.pivot_table(
        bc_sat_df, index=image_meta_col_list + ["Ch"], columns=["type"]
    ).reset_index()
    bc_sat_df.columns = image_meta_col_list + ["Ch", "PercentMax", "StdIntensity"]
    bc_saturation_ymax = max(bc_sat_df.PercentMax)
    if bc_saturation_ymax < 0.2:
        bc_saturation_ymax = 0.2

    for well in bc_sat_df.loc[:, image_cols["well"]].squeeze().unique():
        bc_saturation_gg = (
            gg.ggplot(
                bc_sat_df.loc[bc_sat_df.loc[:, image_cols["well"]] == well, :],
                gg.aes(x="StdIntensity", y="PercentMax", label=image_cols["site"]),
            )
            + gg.coord_fixed(ratio=0.25)
            + gg.geom_text(size=6)
            + gg.facet_wrap("~Ch", ncol=4, scales="free")
            + gg.ylim([0, bc_saturation_ymax])
            + gg.theme_bw()
            + gg.ggtitle(f"Barcoding Image Saturation (Well: {well}) \n {plate}")
            + gg.theme(
                strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
                strip_text=gg.element_text(size=7),
                axis_text=gg.element_text(size=6),
                subplots_adjust={"wspace": 0.7},
            )
        )

        output_file = pathlib.Path(figures_output, f"bc_saturation_{well}.png")
        if check_if_write(output_file, force, throw_warning=True):
            bc_saturation_gg.save(
                output_file,
                dpi=300,
                width=5,
                height=(barcoding_cycles + 2),
                verbose=False,
            )

# Create list of questionable channel correlations (alignments)
correlation_col_prefix = "Correlation_Correlation_"
corr_base_cols = image_meta_col_list.copy() + ["site"]

corr_qc_cols = [col for col in image_df.columns if correlation_col_prefix in col]

if len(corr_qc_cols) > 0:
    image_corr_df = image_df.loc[:, corr_qc_cols + corr_base_cols]

    image_corr_list = []
    for col in corr_qc_cols:
        image_corr_list.append(
            image_corr_df.loc[image_corr_df[col] < correlation_threshold]
        )

    image_corr_df = (
        pd.concat(image_corr_list).drop_duplicates(subset="site").reset_index()
    )

    for col in corr_qc_cols:
        image_corr_df.loc[(image_corr_df[col] >= correlation_threshold), col] = "pass"

    if len(image_corr_df.index) > 0:
        corr_output_file = pathlib.Path(results_output, "flagged_correlations.csv")

        if check_if_write(corr_output_file, force, throw_warning=True):
            image_corr_df.to_csv(corr_output_file)
