import os
import sys
import pathlib
import logging
import traceback
import argparse
import pandas as pd
import plotnine as gg

sys.path.append("config")
from utils import parse_command_args, process_configuration

recipe_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(recipe_path, "scripts"))
from cell_quality_utils import CellQuality
from io_utils import check_if_write, read_csvs_with_chunksize

# Configure logging
logfolder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "logs")
if not os.path.isdir(logfolder):
    os.mkdir(logfolder)
logging.basicConfig(
    filename=os.path.join(logfolder, "4.image-and-segmentation-qc.log"),
    level=logging.INFO,
)


def handle_excepthook(exc_type, exc_value, exc_traceback):
    logging.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))
    traceback_details = "\n".join(traceback.extract_tb(exc_traceback).format())
    print(f"Uncaught Exception: {traceback_details}")


sys.excepthook = handle_excepthook

# Configure experiment
args = parse_command_args()

batch_id = args.batch_id
options_config_file = args.options_config_file
experiment_config_file = args.experiment_config_file
split_step = args.split_step

config = process_configuration(
    batch_id,
    step="preprocess--summarize-plate",
    options_config=options_config_file,
    experiment_config=experiment_config_file,
)
logging.info(f"Config used:{config}")

# Defines the sections of the config file

# Defines the variables set in the config file
barcoding_cycles = config["experiment"]["barcoding_cycles"]
sites_per_image_grid_side = config["experiment"]["sites_per_image_grid_side"]
split_info = config["experiment"]["split"][split_step]

ignore_files = config["options"]["core"]["ignore_files"]
cell_filter = config["options"]["core"]["cell_quality"]["cell_filter"]
quality_func = config["options"]["core"]["cell_quality"]["categorize_cell_quality"]
cell_category_order = config["options"]["core"]["cell_quality"]["cell_category_order"]
cell_category_colors = config["options"]["core"]["cell_quality"]["cell_category_colors"]

input_image_file = config["files"]["image_file"]
cell_count_file = config["files"]["cell_count_file"]
output_resultsdir = config["directories"]["preprocess"]["results"]
output_figuresdir = config["directories"]["preprocess"]["figures"]

spot_config = config["options"]["preprocess"]["process-spots"]
image_cols = spot_config["image_cols"]

plate_summary_config = config["options"]["preprocess"]["summarize-plate"]
correlation_threshold = plate_summary_config["correlation_threshold"]
painting_image_names = plate_summary_config["painting_image_names"]
barcoding_prefix = plate_summary_config["barcoding_prefix"]
force = plate_summary_config["force_overwrite"]
perform = plate_summary_config["perform"]

# check if this step should be performed
if not perform:
    sys.exit("Config file set to perform=False, not performing {}".format(__file__))

# Forced overwrite can be achieved in one of two ways.
# The command line overrides the config file, check here if it is provided
if not force:
    force = args.force

print("Starting 4.image-and-segmentation-qc.")
logging.info(f"Started 4.image-and-segmentation-qc.")

cell_count_df = read_csvs_with_chunksize(cell_count_file, sep="\t")

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
    .assign(site_location=sites_list)
)

# Create total_cell_count
cell_count_bysite_df = (
    cell_count_df.groupby("site")["cell_count"]
    .sum()
    .reset_index()
    .rename(columns={"cell_count": "total_cell_count"})
)

# Add total_cell_count to cell_count_df and add in x, y coordinates for plotting
cell_count_df = cell_count_df.merge(cell_count_bysite_df, on="site").merge(
    loc_df, on="site_location"
)

# Determine list of plates in dataset
platelist = cell_count_df["plate"].unique()

# Plot total number of cells per well
cell_count_totalcells_df = (
    cell_count_df.groupby(["x_loc", "y_loc", "well", "site_location", "site"])[
        "total_cell_count"
    ]
    .mean()
    .reset_index()
)

for plate in platelist:
    os.makedirs(output_figuresdir, exist_ok=True)
    by_well_gg = (
        gg.ggplot(
            cell_count_totalcells_df.loc[
                cell_count_totalcells_df["site"].str.contains(plate)
            ],
            gg.aes(x="x_loc", y="y_loc"),
        )
        + gg.geom_point(gg.aes(fill="total_cell_count"), shape="s", size=6)
        + gg.geom_text(gg.aes(label="site_location"), color="lightgrey", size=6)
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
        + gg.scale_fill_cmap(name="Number of Cells")
    )

    output_file = pathlib.Path(
        output_figuresdir, f"plate_layout_cells_count_per_well_{plate}.png"
    )
    if check_if_write(output_file, force, throw_warning=True):
        by_well_gg.save(output_file, dpi=300, verbose=False)

# Plot cell category ratios per well and empty cells per well
ratio_df = pd.pivot_table(
    cell_count_df,
    values="cell_count",
    index=["site", "plate", "well", "site_location", "x_loc", "y_loc"],
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
)
empty_df = ratio_df.assign(Percent_Empty=ratio_df["Empty"] / ratio_df["Sum"] * 100)
empty_df = (
    empty_df[["Percent_Empty"]]
    .stack()
    .to_frame()
    .reset_index()
    .rename(columns={0: "Percent Empty"})
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
    "Pass_Fail_withempty": "Included/Excluded (with empty cells)",
    "Pass_Fail_0empty": "Included/Excluded (without empty cells)",
}
ratio_df = ratio_df.assign(
    cell_quality_recode=ratio_df.Cell_Quality.replace(quality_recode)
)

for plate in platelist:
    ratio_gg = (
        gg.ggplot(
            ratio_df.loc[ratio_df["plate"].str.contains(plate)],
            gg.aes(x="x_loc", y="y_loc"),
        )
        + gg.geom_point(gg.aes(fill="Ratio"), shape="s", size=6)
        + gg.geom_text(gg.aes(label="site_location"), size=6, color="lightgrey")
        + gg.facet_wrap("~cell_quality_recode + well", ncol=3, scales="fixed")
        + gg.theme_bw()
        + gg.ggtitle(f"Cells Included/Excluded by Cell Quality\n{plate}")
        + gg.theme(
            axis_text=gg.element_blank(),
            axis_title=gg.element_blank(),
            strip_background=gg.element_rect(
                colour="black", fill="#fdfff4", height=1 / 5
            ),
            strip_text=gg.element_text(y=1.1),
        )
        + gg.scale_fill_cmap(name="Ratio")
        + gg.coord_fixed()
    )

    output_file = pathlib.Path(
        output_figuresdir, f"plate_layout_ratios_per_well_{plate}.png"
    )
    if check_if_write(output_file, force, throw_warning=True):
        ratio_gg.save(
            output_file,
            dpi=300,
            width=(len(ratio_df["well"].unique()) + 2),
            height=6,
            verbose=False,
        )

    empty_gg = (
        gg.ggplot(
            empty_df.loc[empty_df["plate"].str.contains(plate)],
            gg.aes(x="x_loc", y="y_loc"),
        )
        + gg.geom_point(gg.aes(fill="Percent Empty"), shape="s", size=6)
        + gg.geom_text(gg.aes(label="site_location"), size=6, color="lightgrey")
        + gg.facet_wrap("~well", ncol=3, scales="fixed")
        + gg.theme_bw()
        + gg.ggtitle(f"Percent Empty Cells\n{plate}")
        + gg.theme(
            axis_text=gg.element_blank(),
            axis_title=gg.element_blank(),
            strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
        )
        + gg.scale_fill_cmap(name="Percent")
        + gg.coord_fixed()
    )

    output_file = pathlib.Path(
        output_figuresdir, f"plate_layout_percent_empty_cells_per_well_{plate}.png"
    )
    if check_if_write(output_file, force, throw_warning=True):
        empty_gg.save(
            output_file, dpi=300, verbose=False,
        )

# Load image file
image_df = read_csvs_with_chunksize(input_image_file, sep="\t")
image_meta_col_list = list(image_cols.values())

# Add in x, y coordinates for plotting
image_df[image_cols["site"]] = image_df[image_cols["site"]].astype(int)
image_df = image_df.merge(
    loc_df, how="left", left_on=image_cols["site"], right_on="site_location"
)

# Plot correlation between BC and CP DAPI images (alignment)
correlation_col_prefix = "Correlation_Correlation_"
for col in image_df.columns:
    if correlation_col_prefix in col:
        if col.count("DAPI") > 1:
            if "Mask" in col:
                corr_qc_col = col
image_df_subset = image_df.dropna(subset=[corr_qc_col])
for plate in platelist:
    correlation_gg = (
        gg.ggplot(
            image_df_subset.loc[
                image_df_subset[image_cols["plate"]].str.contains(plate)
            ],
            gg.aes(x="x_loc", y="y_loc"),
        )
        + gg.geom_point(gg.aes(fill=corr_qc_col), shape="s", size=6)
        + gg.geom_text(gg.aes(label="site_location"), color="lightgrey", size=6)
        + gg.facet_wrap(f"~{image_cols['well']}")
        + gg.coord_fixed()
        + gg.theme_bw()
        + gg.ggtitle(
            f"Correlation between BC & CP DAPI images\n(Alignment quality)\n{plate}"
        )
        + gg.theme(
            axis_text=gg.element_blank(),
            axis_title=gg.element_blank(),
            strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
        )
        + gg.scale_fill_cmap(name="Correlation")
    )

    output_file = pathlib.Path(
        output_figuresdir,
        f"plate_layout_BC_to_CP_DAPI_correlation_per_well_{plate}.png",
    )
    if check_if_write(output_file, force, throw_warning=True):
        correlation_gg.save(output_file, dpi=300, verbose=False)

# Plot final compartment thresholds per well
threshold_col_prefix = "Threshold_FinalThreshold_"
for threshhold_compartment in ["Cells", "Nuclei"]:
    threshold_full_col = f"{threshold_col_prefix}{threshhold_compartment}"

    if threshold_full_col not in image_df.columns:
        continue

    image_df_subset = image_df.dropna(subset=[threshold_full_col])
    for plate in platelist:
        compartment_finalthresh_gg = (
            gg.ggplot(
                image_df_subset.loc[
                    image_df_subset[image_cols["plate"]].str.contains(plate)
                ],
                gg.aes(x="x_loc", y="y_loc"),
            )
            + gg.geom_point(gg.aes(fill=threshold_full_col), shape="s", size=6)
            + gg.geom_text(gg.aes(label="site_location"), color="lightgrey", size=6)
            + gg.facet_wrap(f"~{image_cols['well']}")
            + gg.coord_fixed()
            + gg.theme_bw()
            + gg.ggtitle(f"{threshhold_compartment} Thresholds \n {plate}")
            + gg.theme(
                axis_text=gg.element_blank(),
                axis_title=gg.element_blank(),
                strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
            )
            + gg.scale_fill_cmap(name="Threshold")
            + gg.labs(fill="Threshold")
        )
        output_file = pathlib.Path(
            output_figuresdir,
            f"plate_layout_{threshhold_compartment}_FinalThreshold_per_well_{plate}.png",
        )
        if check_if_write(output_file, force, throw_warning=True):
            compartment_finalthresh_gg.save(output_file, dpi=300, verbose=False)

# Create list of sites with confluent regions
confluent_col = "Math_PercentConfluent"
if confluent_col in image_df.columns:
    # Plot percent confluent regions per well
    image_df_subset = image_df.dropna(subset=[confluent_col])
    for plate in platelist:
        percent_confluent_gg = (
            gg.ggplot(
                image_df_subset.loc[
                    image_df_subset[image_cols["plate"]].str.contains(plate)
                ],
                gg.aes(x="x_loc", y="y_loc"),
            )
            + gg.geom_point(gg.aes(fill=confluent_col), shape="s", size=6)
            + gg.geom_text(gg.aes(label="site_location"), color="lightgrey", size=6)
            + gg.facet_wrap(f"~{image_cols['well']}")
            + gg.coord_fixed()
            + gg.theme_bw()
            + gg.ggtitle(f"Percent of Image That Is Confluent \n {plate}")
            + gg.theme(
                axis_text=gg.element_blank(),
                axis_title=gg.element_blank(),
                strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
            )
            + gg.scale_fill_cmap(name="Percent Confluent")
            + gg.labs(fill="Percent")
        )
        output_file = pathlib.Path(
            output_figuresdir, f"plate_layout_PercentConfluent_per_well_{plate}.png"
        )
        if check_if_write(output_file, force, throw_warning=True):
            percent_confluent_gg.save(output_file, dpi=300, verbose=False)

    confluent_cols = ["site_location", confluent_col] + image_meta_col_list
    confluent_df = image_df.loc[image_df[confluent_col] > 0]
    confluent_df = (
        confluent_df[confluent_cols]
        .sort_values(by=["site_location"])
        .reset_index(drop=True)
    )

    if len(confluent_df.index) > 0:
        confluent_output_file = pathlib.Path(
            output_resultsdir, "sites_with_confluent_regions.csv"
        )
        if check_if_write(confluent_output_file, force, throw_warning=True):
            confluent_df.to_csv(confluent_output_file)

# Power Log Log Slope on Cell Painting images (proxy for focus)
# Any point too high or too low may have focus issues
pll_col_prefix = "ImageQuality_PowerLogLogSlope_"
PLLS_df_cols = image_meta_col_list.copy()
for x in image_df.columns.tolist():
    if pll_col_prefix in x:
        if "Cycle" not in x:
            PLLS_df_cols.append(x)
PLLS_df = image_df.loc[:, PLLS_df_cols]
PLLS_df = PLLS_df.melt(id_vars=image_meta_col_list, var_name="channel").replace(
    {pll_col_prefix: ""}, regex=True
)
for plate in platelist:
    PLLS_gg = (
        gg.ggplot(
            PLLS_df.loc[PLLS_df[image_cols["plate"]].str.contains(plate)],
            gg.aes(x=image_cols["site"], y="value", label=image_cols["site"]),
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

    output_file = pathlib.Path(output_figuresdir, f"PLLS_per_well_{plate}.png")
    if check_if_write(output_file, force, throw_warning=True):
        PLLS_gg.save(
            output_file,
            width=(len(PLLS_df[image_cols["well"]].unique()) + 2),
            height=8,
            dpi=300,
            verbose=False,
        )

# Outputs list of sites that are saturated in any channel
# Cell Painting images use >1% saturated, Barcoding images uses >.2% saturated
saturated_col_prefix = "ImageQuality_PercentMaximal_"
cp_sat_cols = []
bc_sat_cols = []
nts = ["A", "C", "G", "T"]
cp_sat_df = pd.DataFrame()
bc_sat_df = pd.DataFrame()
sat_df = pd.DataFrame()

# Create lists of columns measuring saturation
for col in image_df.columns:
    if saturated_col_prefix in col:
        if "Cycle" not in col:
            cp_sat_cols.append(col)
        if "Cycle" in col:
            bc_sat_cols.append(col)

# Create df of sites that fail CP saturation
if all(x in image_df.columns.tolist() for x in cp_sat_cols):
    for col in cp_sat_cols:
        temp = image_df[image_df[col] > 1]
        cp_sat_df = cp_sat_df.append(temp)
    if len(cp_sat_df) != 0:
        cp_sat_df.loc[:, "Fails_CP_Sat"] = "Fails"

# Create df of sites that fail BC saturation
if all(x in image_df.columns.tolist() for x in bc_sat_cols):
    for col in bc_sat_cols:
        temp = image_df[image_df[col] > 0.2]
        bc_sat_df = bc_sat_df.append(temp)
    if len(bc_sat_df) != 0:
        bc_sat_df.loc[:, "Fails_BC_Sat"] = "Fails"

if len(cp_sat_df) > 0 and len(bc_sat_df) > 0:
    sat_df = (
        cp_sat_df.set_index("Metadata_Site_Full")
        .combine_first(bc_sat_df.set_index("Metadata_Site_Full"))
        .reset_index()
    )
    addn_cols = ["Metadata_Site_Full", "Fails_CP_Sat", "Fails_BC_Sat"]
elif len(cp_sat_df) > 0:
    sat_df = cp_sat_df
    addn_cols = ["Metadata_Site_Full", "Fails_CP_Sat"]
elif len(bc_sat_df) > 0:
    sat_df = bc_sat_df
    addn_cols = ["Metadata_Site_Full", "Fails_BC_Sat"]

if not sat_df.empty:
    sat_df_cols = cp_sat_cols + bc_sat_cols + addn_cols
    sat_df = sat_df.loc[:, sat_df_cols]
    if len(cp_sat_df) > 0:
        sat_df["Fails_CP_Sat"].fillna("Passes", inplace=True)
    if len(bc_sat_df) > 0:
        sat_df["Fails_BC_Sat"].fillna("Passes", inplace=True)

# saturated_sites.csv does not save if empty
if len(sat_df.index) > 0:
    sat_output_file = pathlib.Path(output_resultsdir, "saturated_sites.csv")
    if check_if_write(output_file, force, throw_warning=True):
        sat_df.to_csv(sat_output_file)

# Plots saturation in Cell Painting images
# x = std dev of intensity (to find images that have unusually bright spots)
# y = % image that is saturated (to find images that are unusually bright)
# Look at points off cluster where y > 1
intensity_col_prefix = "ImageQuality_StdIntensity_"
cp_sat_df_cols = image_meta_col_list.copy()
for x in image_df.columns.tolist():
    if "Cycle" not in x:
        if intensity_col_prefix in x:
            cp_sat_df_cols.append(x)
        if saturated_col_prefix in x:
            cp_sat_df_cols.append(x)

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
    for plate in platelist:
        cp_saturation_gg = (
            gg.ggplot(
                cp_sat_df.loc[cp_sat_df[image_cols["plate"]].str.contains(plate)],
                gg.aes(x="StdIntensity", y="PercentMax", label=image_cols["site"]),
            )
            + gg.geom_text(size=6)
            + gg.facet_wrap(
                ["Ch", image_cols["well"]],
                nrow=len(painting_image_names),
                scales="fixed",
            )
            + gg.ylim([0, cp_saturation_ymax])
            + gg.theme_bw()
            + gg.ggtitle(f"Cell Painting Image Saturation \n {plate}")
            + gg.theme(
                strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
                strip_text=gg.element_text(size=7),
                axis_text=gg.element_text(size=6),
                axis_text_x=gg.element_text(rotation=90),
                subplots_adjust={"wspace": 0.2},
            )
        )
        output_file = pathlib.Path(output_figuresdir, f"cp_saturation_{plate}.png")
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
for x in image_df.columns.tolist():
    if "Cycle" in x:
        if intensity_col_prefix in x:
            bc_sat_df_cols.append(x)
        if saturated_col_prefix in x:
            bc_sat_df_cols.append(x)

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

    ch_order = []
    for ch in bc_sat_df["Ch"].unique():
        if "DAPI" not in ch:
            ch_order.append(ch)
    for ch in bc_sat_df["Ch"].unique():
        if "DAPI" in ch:
            ch_order.append(ch)
    bc_sat_df["Ch"] = bc_sat_df.Ch.astype("category")
    bc_sat_df["Ch"] = bc_sat_df["Ch"].cat.reorder_categories(ch_order)

    bc_sat_df.columns = image_meta_col_list + ["Ch", "PercentMax", "StdIntensity"]
    bc_saturation_ymax = max(bc_sat_df.PercentMax)
    if bc_saturation_ymax < 0.2:
        bc_saturation_ymax = 0.2

    for well in bc_sat_df.loc[:, image_cols["well"]].squeeze().unique():
        for plate in platelist:
            bc_saturation_gg = (
                gg.ggplot(
                    bc_sat_df.loc[
                        (bc_sat_df[image_cols["well"]] == well)
                        & (bc_sat_df[image_cols["plate"]] == plate)
                    ],
                    gg.aes(x="StdIntensity", y="PercentMax", label=image_cols["site"]),
                )
                + gg.geom_text(size=6)
                + gg.facet_wrap("~Ch", ncol=4, scales="fixed")
                + gg.ylim([0, bc_saturation_ymax])
                + gg.theme_bw()
                + gg.ggtitle(f"Barcoding Image Saturation (Well: {well}) \n {plate}")
                + gg.theme(
                    strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
                    strip_text=gg.element_text(size=7),
                    axis_text=gg.element_text(size=6),
                    axis_text_x=gg.element_text(rotation=90),
                    subplots_adjust={"wspace": 0.2},
                )
            )

            output_file = pathlib.Path(
                output_figuresdir, f"bc_saturation_{well}_{plate}.png"
            )
            if check_if_write(output_file, force, throw_warning=True):
                bc_saturation_gg.save(
                    output_file,
                    dpi=300,
                    width=5,
                    height=(barcoding_cycles + 2),
                    verbose=False,
                )
print("Finished 4.image-and-segmentation-qc.")
logging.info(f"Finished 4.image-and-segmentation-qc.")
