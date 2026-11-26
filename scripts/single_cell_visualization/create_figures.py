"""
Make single-cell montage figures for selected genes.

Usage:
    python create_figures.py --data-set NAME --genes GENE [GENE ...] \
        --grouping {gene,guide} --method {random,representative} \
        [--defaults-config DEFAULTS_CONFIG] [--experiment-config EXPERIMENT_CONFIG] \
        [--n-cells N] [--box-size PX] [--dpi DPI] [--percentile P] \
        [--no-objects] [--no-color] [--no-bc-nuclei] [--bc-nuclei-name NAME]

Required: --data-set, --genes, --grouping, --method.
Optional: --defaults-config and --experiment-config (default to config/defaults.yaml and
config/experiment.json, relative to the working directory, as in run.py) and figure
parameters (defaults shown in --help).

Reads:
- {data}/Visualization_metadata.json (written by 0.image-qc.py)
- {single_cell}/{plate}_alignments_{data_set_name}.csv (written by 0.image-qc.py)
- {single_cell}/single_cell_by_gene/{data_set_name}_{grouping}_normalized_{gene}.csv.gz (written by 5.normalize.py)
- {plate}_load_data_pipeline9.csv (see LOAD_DATA_TEMPLATE below)

Outputs the following files:
- outputs/{gene}_{data_set_name}_{cell_selection_method}_{plate}.png

ASSUMES:
- FileName_<channel> columns of load_data hold the full path to each image
- Objects are stored in structure OBJECT_ROOT/BATCH/PLATE-WELL-SITE/segmentation_masks/Plate_PLATE_Well_WELL_Site_SITE_Cells_Objects.tiff

BEFORE RUNNING:
- Run the recipe through 5.normalize.py to generate the single-cell CSVs
- Download load_data.csv's used for analysis to load_data/
"""

import argparse
import os
import sys
import json
import glob
import pandas as pd
import matplotlib.pyplot as plt

this_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(this_dir)  # so `from sc_vis_utils import ...` works from any cwd
recipe_path = os.path.dirname(os.path.dirname(this_dir))
sys.path.append(os.path.join(recipe_path, "utils"))

from sc_vis_utils import filter_out_edge_single_cells, extract_single_cell_samples, make_figure
from io_utils import load_configs, parse_data_set

#################################
# DEFAULT FIGURE PARAMETERS (override on the command line)
#################################

DEFAULT_N_CELLS = 20  # number of cells to return
DEFAULT_BOX_SIZE = 100  # in pixels, size of bounding box to crop around cells
DEFAULT_DPI = 300
DEFAULT_BC_NUCLEI_NAME = "Cycle01_DAPI"
DEFAULT_PERCENTILE = 98  # percentile for brightness adjustment
DEFAULT_OBJECT_ROOT = None  # root directory for cell objects


# Channel -> display color, as understood by make_figure.CP_to_RGB_single.
# Edit if your phenotyping channels are not listed.
CHANNEL_COLORS = {
    "DNA": "Blue",
    "ER": "Green",
    "RNA": "Yellow",
    "AGP": "Red",
    "Mito": "Magenta",
    "DAPI": "Blue",
    "WGA": "Red",
    "Phalloidin": "Green",
    "Brightfield": "White",
}

outdir = "single_cell_figures"  # location to output figures


#################################
# FUNCTIONS
#################################


def process_figures(
    path_to_defaults_config,
    path_to_experiment_config,
    data_set,
    genelist,
    grouping,
    cell_selection_method,
    n_cells=DEFAULT_N_CELLS,
    box_size=DEFAULT_BOX_SIZE,
    dpi_quality=DEFAULT_DPI,
    objects=True,  # add cell objects
    color=True,  # add color composite
    bc_nuclei=True,  # add nuclei from SBS images
    bc_nuclei_name=DEFAULT_BC_NUCLEI_NAME,
    percentile=DEFAULT_PERCENTILE,
    object_root=DEFAULT_OBJECT_ROOT,
):
    defaults_config, experiment_config = load_configs(
        path_to_defaults_config, path_to_experiment_config
    )
    data_sets = experiment_config["data_sets"]
    if data_set not in data_sets:
        raise ValueError(
            f"--data-set {data_set!r} not found in experiment config. Available: {list(data_sets.keys())}"
        )

    out_root = defaults_config["directory_structure"]["root"]
    data_folder = defaults_config["directory_structure"]["data"]
    single_cell = defaults_config["directory_structure"]["single_cell"]
    single_cell_dir = os.path.join(out_root, single_cell)
    sc_files_dir = os.path.join(single_cell_dir, "single_cell_by_gene")
    file_naming = "{data_set_name}_{grouping}_normalized_{gene}.csv.gz"

    with open(
        os.path.join(out_root, data_folder, "Visualization_metadata.json"),
        "r",
        encoding="utf-8",
    ) as f:
        vis_metadata = json.load(f)
    im_height = vis_metadata["image_height"]
    im_width = vis_metadata["image_width"]
    chans = vis_metadata["inferred_phenotyping_channels"]
    im_dict = {c: CHANNEL_COLORS.get(c, "White") for c in chans}

    if objects and not object_root:
        raise ValueError("--object-root is required unless --no-objects is passed.")

    os.makedirs(outdir, exist_ok=True)

    # load all load_data from {recipe}/load_data and combine
    load_data_files = sorted(
        glob.glob(os.path.join(recipe_path, "load_data", "**", "*.csv"), recursive=True)
    )
    if not load_data_files:
        raise FileNotFoundError(
            f"No load_data csvs found in {os.path.join(recipe_path, 'load_data')}"
        )
    all_loaddata = pd.concat(
        [pd.read_csv(file) for file in load_data_files], ignore_index=True
    )

    for data_set_name in [data_set]:
        print(f"Processing treatment: {data_set_name}")
        batch, plate_list, list_plate_well_tuples = parse_data_set(
            data_sets[data_set_name]
        )
        for plate in plate_list:
            # Read in metadata output by recipe
            alignments = pd.read_csv(
                os.path.join(single_cell_dir, f"{plate}_alignments_{data_set_name}.csv")
            )
            alignments["Metadata_Site"] = alignments["Metadata_Site"].astype(int)
            # The alignments csv always exists, but only has alignment columns
            # (containing _X / _Y) when image alignment shifts were measured.
            align_cols = [c for c in alignments.columns if not c.startswith("Metadata")]
            use_alignments = any("_X" in c for c in align_cols) and any(
                "_Y" in c for c in align_cols
            )
            if not use_alignments:
                print(f"No alignment columns for {plate} in {data_set_name}; not applying offsets.")

            # this plate's load_data, with path to objects output during analysis
            loaddata = all_loaddata[all_loaddata["Metadata_Plate"] == plate].copy()
            loaddata["Path_CellObjects"] = (
                f"{object_root}/{batch}/{plate}-"
                + loaddata["Metadata_Well"]
                + "-"
                + loaddata["Metadata_Site"].astype(str)
                + f"/segmentation_masks/Plate_{plate}_Well_"
                + loaddata["Metadata_Well"]
                + "_Site_"
                + loaddata["Metadata_Site"].astype(str)
                + "_Cells_Objects.tiff"
            )

            for gene in genelist:
                print(f"Making figure for gene/guide {gene}")
                # load the single cell profiles for the specified gene and data set
                df = pd.read_csv(
                    os.path.join(
                        sc_files_dir,
                        file_naming.format(
                            data_set_name=data_set_name,
                            grouping=grouping,
                            gene=gene,
                        ),
                    )
                )
                # clean up the dataframe
                df["Metadata_Site"] = df["Metadata_Site"].astype(int)
                df["Metadata_Label"] = (
                    df["Metadata_Well"] + "-" + df["Metadata_Site"].astype(str)
                )

                # add the alignment columns, only if there are any
                if use_alignments:
                    df = df.merge(
                        alignments,
                        on=["Metadata_Plate", "Metadata_Well", "Metadata_Site"],
                        how="left",
                    )
                df = df.merge(
                    loaddata,
                    on=["Metadata_Plate", "Metadata_Well", "Metadata_Site"],
                    how="left",
                )

                # filter out cells that are too close to the edge of the image
                df, _ = filter_out_edge_single_cells.edgeCellFilter(
                    df,
                    "Nuclei_AreaShape_Center_X",
                    "Nuclei_AreaShape_Center_Y",
                    im_height,
                    im_width,
                    box_size / 2,
                )

                if df.shape[0] > 0:
                    (
                        df,
                        cp_features_analysis,
                    ) = extract_single_cell_samples.extract_single_cell_samples(
                        df,
                        n_cells,
                        cell_selection_method,
                        platefilter=[plate],
                    )
                    print(
                        f"Filtered to {len(df)} cells using {cell_selection_method} selection method."
                    )
                else:
                    print(f"No cells found for {gene} after filtering.")
                    continue

                figure = make_figure.make_figure(
                    gene,
                    df,
                    im_dict,
                    cell_selection_method,
                    box_size,
                    objects,
                    color,
                    bc_nuclei,
                    bc_nuclei_name,
                    percentile,
                    use_alignments=use_alignments,
                )

                figure.savefig(
                    os.path.join(
                        outdir,
                        f"{gene}_{data_set_name}_{cell_selection_method}_{plate}.png",
                    ),
                    dpi=dpi_quality,
                )
                plt.close("all")


#################################
# MAIN USER INTERACTION
#################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.split("Usage:")[0].strip())
    parser.add_argument(
        "--defaults-config",
        dest="defaults_config",
        default="config/defaults.yaml",
        help="Path to defaults config file (default: %(default)s)",
    )
    parser.add_argument(
        "--experiment-config",
        dest="experiment_config",
        default="config/experiment.json",
        help="Path to experiment config file (default: %(default)s)",
    )

    req = parser.add_argument_group("required parameters")
    req.add_argument("--data-set", required=True, help="Data set name; must match a name in the experiment config")
    req.add_argument("--genes", required=True, nargs="+", help="Gene(s) (or guide(s)) to make figures for")
    req.add_argument("--grouping", required=True, choices=["gene", "guide"])
    req.add_argument("--method", required=True, choices=["random", "representative"])

    opt = parser.add_argument_group("optional figure parameters")
    opt.add_argument("--n-cells", type=int, default=DEFAULT_N_CELLS, help="Number of cells to return (default: %(default)s)")
    opt.add_argument("--box-size", type=int, default=DEFAULT_BOX_SIZE, help="Pixel size of bounding box cropped around cells (default: %(default)s)")
    opt.add_argument("--dpi", type=int, default=DEFAULT_DPI, help="Figure dpi (default: %(default)s)")
    opt.add_argument("--percentile", type=float, default=DEFAULT_PERCENTILE, help="Percentile for brightness adjustment (default: %(default)s)")
    opt.add_argument("--bc-nuclei-name", default=DEFAULT_BC_NUCLEI_NAME, help="SBS nuclei image name (default: %(default)s)")
    opt.add_argument("--object-root", default=DEFAULT_OBJECT_ROOT, help="Root directory for cell objects (default: %(default)s)")
    opt.add_argument("--no-objects", dest="objects", action="store_false", help="Omit cell objects")
    opt.add_argument("--no-color", dest="color", action="store_false", help="Omit color composite")
    opt.add_argument("--no-bc-nuclei", dest="bc_nuclei", action="store_false", help="Omit SBS nuclei")

    args = parser.parse_args()
    process_figures(
        args.defaults_config,
        args.experiment_config,
        data_set=args.data_set,
        genelist=args.genes,
        grouping=args.grouping,
        cell_selection_method=args.method,
        n_cells=args.n_cells,
        box_size=args.box_size,
        dpi_quality=args.dpi,
        objects=args.objects,
        color=args.color,
        bc_nuclei=args.bc_nuclei,
        bc_nuclei_name=args.bc_nuclei_name,
        percentile=args.percentile,
        object_root=args.object_root,
    )
