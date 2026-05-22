"""
Does not allow for selection of overwrite.
"""

import os
import sys
import logging
import traceback
import math
import pandas as pd
import seaborn.objects as so
import matplotlib.pyplot as plt


recipe_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(recipe_path, "utils"))
from io_utils import load_configs, parse_data_set, read_csvs_with_chunksize
from plotting_utils import make_loc_df, make_plate_layout_plots

# Configure logging
logfolder = os.path.join(recipe_path, "logs")
if not os.path.isdir(logfolder):
    os.mkdir(logfolder)
logging.basicConfig(
    filename=os.path.join(logfolder, "0.image-qc.log"),
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
            f"Allowed skip limit of {allowed_skips} reached. Stopping 0.image-qc.",
            type="warning",
        )
        sys.exit(1)
    else:
        return allowed_skip_counter


def make_intensity_plot(images_df, cols, outpath):
    df = images_df.set_index(["Metadata_Identifier", "Metadata_Site", "Metadata_Well"])[
        cols
    ]
    df.columns = df.columns.str.replace("ImageQuality_", "")
    df.columns = df.columns.str.split("_", n=1, expand=True)
    df = df.stack(level=1, future_stack=True)
    df = df.reset_index().rename(columns={"level_3": "Channel"})
    p = (
        so.Plot(df, x="PercentMaximal", y="StdIntensity", text="Metadata_Site")
        .facet(col="Metadata_Well", row="Channel")
        .add(so.Text())
        .label(
            y="Percent Image Saturated",
            x="StdDev of Intensity (unusally bright spots)",
        )
    ).theme({"axes.facecolor": "w", "axes.edgecolor": "black"})
    p.save(outpath, dpi=300)


def process_qc(path_to_defaults_config, path_to_experiment_config):
    printandlog("Started 0.image-qc.")
    defaults_config, experiment_config = load_configs(
        path_to_defaults_config, path_to_experiment_config
    )

    # define experiment variables
    file_location = experiment_config["file_location"]
    data_sets = experiment_config["data_sets"]
    compartments = experiment_config["compartments"]

    # define defaults variables
    out_root = defaults_config["directory_structure"]["root"]
    data_folder = defaults_config["directory_structure"]["data"]
    single_cell = defaults_config["directory_structure"]["single_cell"]
    figures_folder = defaults_config["directory_structure"]["figures"]
    allowed_skips = defaults_config["process"]["qc"]["allowed_skips"]
    infer_empty_sites = defaults_config["process"]["qc"]["infer_empty_sites"]
    stack_alignment_chans = defaults_config["process"]["qc"]["stack_alignment_chans"]

    outdir_figs = os.path.join(out_root, figures_folder)
    if not os.path.exists(outdir_figs):
        os.makedirs(outdir_figs, exist_ok=True)
    outdir_data = os.path.join(out_root, data_folder)
    if not os.path.exists(outdir_data):
        os.makedirs(outdir_data, exist_ok=True)
    outdir_single_cell = os.path.join(out_root, single_cell)
    if not os.path.exists(outdir_single_cell):
        os.makedirs(outdir_single_cell, exist_ok=True)

    allowed_skip_counter = 0
    for data_set_name in data_sets.keys():
        printandlog(f"Starting {data_set_name}")
        batch, plate_list, list_plate_well_tuples = parse_data_set(
            data_sets[data_set_name]
        )
        try:
            folderlist = os.listdir(os.path.join(file_location, batch))
        except:
            printandlog(f"Failed to find list of folders in {os.path.join(file_location, batch)}. Likely need to correct 'file_location' in config.")
        inferred_empty_sites = []
        all_image_dfs = []
        for plate in plate_list:
            all_align_dfs = []
            for plate_well_tuple in [x for x in list_plate_well_tuples if plate in x]:
                well = plate_well_tuple[1]
                usefolders = [x for x in folderlist if f"{plate}-{well}" in x]
                printandlog(f"Now processing {data_set_name}: {batch}, {plate}, {well}")
                printandlog(f"{len(usefolders)} site folders to process")

                for plate_well_site_folder in usefolders:
                    site = plate_well_site_folder.rsplit("-", 1)[1]
                    try:
                        image_file = os.path.join(
                            file_location,
                            batch,
                            plate_well_site_folder,
                            "Image.csv",
                        )
                        image_df = read_csvs_with_chunksize(image_file)
                        keep_cols = (
                            ["Count_Cells", "Math_PercentConfluent"]
                            + [
                                x
                                for x in image_df.columns
                                if stack_alignment_chans[0] in x
                                and stack_alignment_chans[1] in x
                                and "Correlation_Correlation_" in x
                            ]
                            + [
                                x
                                for x in image_df.columns
                                if "Threshold_FinalThreshold_" in x
                            ]
                            + [
                                x
                                for x in image_df.columns
                                if "ImageQuality_PowerLogLogSlope_" in x
                                and "Cycle" not in x
                            ]
                            + [
                                x
                                for x in image_df.columns
                                if "ImageQuality_PercentMaximal_" in x
                            ]
                            + [
                                x
                                for x in image_df.columns
                                if "ImageQuality_StdIntensity_" in x
                            ]
                            + [x for x in image_df.columns if "Align_" in x]
                        )
                        image_df = (
                            image_df[keep_cols]
                            .assign(
                                Metadata_Batch=batch,
                                Metadata_Plate=plate,
                                Metadata_Well=well,
                                Metadata_Site=site,
                                Metadata_Identifier=f"{batch}-{plate}-{well}-{site}",
                                Metadata_Dataset_Split=data_set_name,
                            )
                            .astype({"Metadata_Site": int})
                        )

                        # Collate all image alignment for by-gene output
                        all_align_dfs.append(
                            image_df[
                                [x for x in image_df.columns if "Align_" in x]
                                + [x for x in image_df.columns if "Metadata" in x]
                            ]
                        )
                        # Collate all per-site Image.csv's for data visualization
                        all_image_dfs.append(image_df)
                    except:
                        if infer_empty_sites:
                            compartment_file_list = [
                                image_file.replace("Image", compartment)
                                for compartment in compartments
                            ]
                            if any(os.path.exists(f) for f in compartment_file_list):
                                allowed_skip_counter = fail_site(
                                    plate_well_site_folder,
                                    image_file,
                                    allowed_skip_counter,
                                    allowed_skips,
                                )
                                continue
                            else:
                                printandlog(
                                    f"Inferred empty site at {plate_well_site_folder}. Continuing without counting as error."
                                )
                                inferred_empty_sites.append(plate_well_site_folder)
                                continue
                        else:
                            allowed_skip_counter = fail_site(
                                plate_well_site_folder,
                                image_file,
                                allowed_skip_counter,
                                allowed_skips,
                            )
                            continue

            # Save Alignment_* columns by plate for by-gene output of normalized profiles
            outfile = os.path.join(
                outdir_single_cell, f"{plate}_alignments_{data_set_name}.csv"
            )
            pd.concat(all_align_dfs).to_csv(outfile, index=False)

        # Add in x, y coordinates to images_df for plotting
        images_df = pd.concat(all_image_dfs)

        sites_per_image_grid_side = int(
            math.ceil(math.sqrt(images_df["Metadata_Site"].astype(int).max()))
        )
        loc_df = make_loc_df(sites_per_image_grid_side)

        images_df = images_df.merge(loc_df, on="Metadata_Site")

        # Plot Cell Count plate layout visualization
        try:
            outpath = os.path.join(
                outdir_figs, f"plate_layout_cells_count_{data_set_name}.png"
            )
            make_plate_layout_plots(images_df, "Count_Cells", "Cell Count", outpath)
        except:
            printandlog(f"Failed to create Cell Count plot for {data_set_name}")

        # Plot confluent regions plate layout visualization
        try:
            outpath = os.path.join(
                outdir_figs, f"plate_layout_image_stack_alignment_{data_set_name}.png"
            )
            make_plate_layout_plots(
                images_df, "Math_PercentConfluent", "Percent Confluent Regions", outpath
            )
        except:
            printandlog(f"Failed to create confluent regions plot for {data_set_name}")

        # Plot image stack alignment plate layout visualization
        try:
            outpath = os.path.join(
                outdir_figs, f"plate_layout_image_stack_alignment_{data_set_name}.png"
            )
            stack_alignment_col = [
                x for x in image_df.columns if "Correlation_Correlation_" in x
            ][0]
            make_plate_layout_plots(
                images_df,
                stack_alignment_col,
                "Correlation between SBS and Phenotyping",
                outpath,
                legend=True,
            )
        except:
            printandlog(
                f"Failed to create image stack alignment plot for {data_set_name}"
            )

        # Plot compartment thresholds plate layout visualization
        for compartment in compartments:
            try:
                outpath = os.path.join(
                    outdir_figs,
                    f"plate_layout_compartment_threshold_{compartment}_{data_set_name}.png",
                )
                make_plate_layout_plots(
                    images_df,
                    f"Threshold_FinalThreshold_{compartment}",
                    f"{compartment} threshold",
                    outpath,
                    legend=True,
                )
            except:
                printandlog(
                    f"Failed to create compartment threshold plot for {compartment} in {data_set_name}"
                )
                printandlog(
                    f"Failure is expected for compartments defined outside of IdentifyPrimaryObjects and IdentifySecondaryObjects modules.\n Cytoplasm is one such example as it is typically a tertiary object.",
                )

        # Power Log Log Slope on Cell Painting images (proxy for focus)
        # Any point too high or too low may have focus issues
        try:
            df = images_df.set_index(
                ["Metadata_Identifier","Metadata_Plate","Metadata_Site", "Metadata_Well"]
            )[[x for x in image_df.columns if "ImageQuality_PowerLogLogSlope_" in x]]
            df.columns = df.columns.str.replace("ImageQuality_", "")
            df.columns = df.columns.str.split("_", n=1, expand=True)
            df = df.stack(level=1, future_stack=True)
            df = df.reset_index().rename(columns={"level_4": "Channel"})
            df['Plate|Well'] = df['Metadata_Plate'] + "|" + df['Metadata_Well'].astype(str)
            outpath = os.path.join(outdir_figs, f"Image_Focus_{data_set_name}.png")
            p = (
                (
                    so.Plot(
                        df,
                        x="Metadata_Site",
                        y="PowerLogLogSlope",
                        text="Metadata_Site",
                    )
                    .facet(row=("Plate|Well"), col="Channel")
                    .add(so.Text())
                )
                .theme({"axes.facecolor": "w", "axes.edgecolor": "black"})
                .plot()
            )
            p.save(outpath, dpi=300)
        except:
            printandlog(f"Failed to create Power Log Log Slope plot in {data_set_name}")

        # Create list of sites with confluent regions
        try:
            confluent_df = image_df.loc[image_df["Math_PercentConfluent"] > 0][
                ["Math_PercentConfluent", "Metadata_Identifier"]
            ].sort_values(by=["Metadata_Identifier"])
            if len(confluent_df) > 0:
                confluent_df.to_csv(
                    os.path.join(
                        out_root,
                        data_folder,
                        f"Sites_With_Confluent_Regions_{data_set_name}.csv",
                    ),
                    index=False,
                )
        except:
            printandlog(f"Failed to create confluent regions .csv for {data_set_name}")

        # Outputs list of sites that are saturated in any channel
        # Cell Painting images use >1% saturated, Barcoding images uses >.2% saturated
        try:
            sat_SBS_df = image_df[
                [
                    x
                    for x in image_df.columns
                    if "ImageQuality_PercentMaximal_" in x and "Cycle" in x
                ]
                + ["Metadata_Identifier"]
            ].set_index("Metadata_Identifier")
            SBS_sat_list = (
                sat_SBS_df.loc[:, (sat_SBS_df > 0.002).any()]
                .reset_index()["Metadata_Identifier"]
                .tolist()
            )
            sat_SBS_df = pd.DataFrame(
                columns=["Metadata_Identifier", "SBS_saturated"],
                data=list(zip(SBS_sat_list, ["Fails"] * len(SBS_sat_list))),
            )

            sat_phenotype_df = image_df[
                [
                    x
                    for x in image_df.columns
                    if "ImageQuality_PercentMaximal_" in x and "Cycle" not in x
                ]
                + ["Metadata_Identifier"]
            ].set_index("Metadata_Identifier")
            phenotype_sat_list = (
                sat_phenotype_df.loc[:, (sat_phenotype_df > 0.01).any()]
                .reset_index()["Metadata_Identifier"]
                .tolist()
            )
            sat_phenotype_df = pd.DataFrame(
                columns=["Metadata_Identifier", "Phenotyping_saturated"],
                data=list(zip(phenotype_sat_list, ["Fails"] * len(phenotype_sat_list))),
            )

            sat_df = sat_SBS_df.merge(
                sat_phenotype_df, on="Metadata_Identifier", how="outer"
            )
            sat_df["Phenotyping_saturated"] = sat_df["Phenotyping_saturated"].fillna(
                "Passes"
            )
            sat_df["SBS_saturated"] = sat_df["SBS_saturated"].fillna("Passes")

            if len(sat_df.index) > 0:
                sat_df.to_csv(
                    os.path.join(
                        out_root,
                        data_folder,
                        f"Sites_With_Saturation_{data_set_name}.csv",
                    ),
                    index=False,
                )
        except:
            printandlog(f"Failed to create saturated sites .csv for {data_set_name}")

        # Plot image intensity plots
        try:
            if len(sat_phenotype_df) > 0:
                try:
                    cols = [
                        x
                        for x in image_df.columns
                        if "ImageQuality_StdIntensity_" in x and "Cycle" not in x
                    ] + [
                        x
                        for x in image_df.columns
                        if "ImageQuality_PercentMaximal_" in x and "Cycle" not in x
                    ]
                    outpath = os.path.join(
                        outdir_figs, f"Image_Saturation_phenotyping_{data_set_name}.png"
                    )
                    make_intensity_plot(images_df, cols, outpath)
                except:
                    printandlog(
                        f"Failed to create phenotyping image intensity plot in {data_set_name}"
                    )
        except:
            continue
        try:
            if len(sat_SBS_df) > 0:
                try:
                    cycles = [
                        x.split("Cycle")[1].split("_")[0]
                        for x in image_df.columns
                        if "ImageQuality_StdIntensity_Cycle" in x
                    ]
                    for cycle in cycles:
                        cols = [
                            x
                            for x in image_df.columns
                            if "ImageQuality_StdIntensity_" in x
                            and f"Cycle{cycle}" in x
                        ] + [
                            x
                            for x in image_df.columns
                            if "ImageQuality_PercentMaximal_" in x
                            and f"Cycle{cycle}" in x
                        ]
                        outpath = os.path.join(
                            outdir_figs,
                            f"Image_Saturation_SBS_Cycle{cycle}_{data_set_name}.png",
                        )
                        make_intensity_plot(images_df, cols, outpath)
                except:
                    printandlog(
                        f"Failed to create SBS image intensity plots in {data_set_name}"
                    )
        except:
            continue

        if inferred_empty_sites:
            inferred_empty_sites_df = pd.DataFrame(
                inferred_empty_sites, columns=["Inferred_Empty_Site_Folder"]
            )
            inferred_empty_sites_df.to_csv(
                os.path.join(
                    out_root,
                    data_folder,
                    f"Inferred_Empty_Sites_{data_set_name}.csv",
                ),
                index=False,
            )

    print("Finished 0.image-qc.")


#################################
# MAIN USER INTERACTION
#################################

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Pass paths to defaults config and experiment config.")
        sys.exit()
    path_to_defaults_config = sys.argv[1]
    path_to_experiment_config = sys.argv[2]
    process_qc(path_to_defaults_config, path_to_experiment_config)
