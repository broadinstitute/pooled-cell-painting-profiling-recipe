# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "marimo>=0.25.2",
#     "matplotlib==3.11.2",
#     "numpy==2.5.4",
#     "obstore==0.11.1",
#     "pandas==3.1.1",
# ]
# ///

import marimo

__generated_with = "0.25.2"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _():
    import os
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np
    import json
    import glob

    return glob, json, os, pd


@app.cell
def _(json):
    public_root = 'scripts/public'
    exp_config = json.load(open('scripts/public/data_sets.json'))
    experiment_list = list(exp_config.keys())
    return exp_config, experiment_list, public_root


@app.cell
def _(mo):
    mo.md("""
    # OPS Experimental Results Explorer v0.1

    ## Broad Institute Imaging Platform, 2026

    ## Start by choosing the sub-experiment
    """)
    return


@app.cell
def _(experiment_list, mo):
    exp_el = mo.ui.dropdown(options = experiment_list, value = experiment_list[0], label = "Choose the sub-experiment to view. Note not all outputs are subdivided by experiment.")
    exp_el
    return (exp_el,)


@app.cell
def _(exp_config, exp_el, mo):
    exp = exp_el.value
    batch_list = list(exp_config[exp].keys())
    if len(batch_list) >1:
        batch_el = mo.ui.dropdown(options = batch_list, value = batch_list[0], label = "Choose the batch of the experiment. Note not all outputs are subdivided by batch.")
        batch_el
    return batch_el, batch_list, exp


@app.cell
def _(batch_el, batch_list, exp, exp_config, mo):
    if len(batch_list) >1:
        batch = batch_el.value
    else:
        batch = batch_list[0]

    plate_list = list(exp_config[exp][batch].keys())
    if len(plate_list) >1:
        plate_el = mo.ui.dropdown(options = plate_list, value = plate_list[0], label = "Choose the batch of the experiment. Note not all outputs are subdivided by plates.")
        plate_el
    return batch, plate_el, plate_list


@app.cell
def _(batch, exp, exp_config, plate_el, plate_list):
    if len(plate_list) >1:
        plate = plate_el.value
    else:
        plate = plate_list[0]
    well_list = exp_config[exp][batch][plate]
    return plate, well_list


@app.cell
def _(mo):
    mo.md("""
    ## Profiling stats
    """)
    return


@app.cell
def _(exp, mo, os, plate, public_root):
    mo.vstack([
        mo.md("### mAP plot"),
        mo.image(os.path.join(public_root,f"data_exploration_{plate}_{exp}_guide_normalized_all_feature_selected_mAP_activity_scatterplot.png"), width = 600)
    ]
    )
    return


@app.cell
def _(exp, mo, os, pd, plate, public_root):
    mAP_df = pd.read_csv(os.path.join(public_root,f"data_exploration_{plate}_{exp}_guide_normalized_all_feature_selected_mAP_activity.csv.gz"))
    mo.vstack([
        mo.md("### mAP Dataframe"),
        mAP_df
    ]
    )
    return


@app.cell
def _(exp, mo, os, plate, public_root):
    mo.vstack([
        mo.md("### Heatmap - per gene"),
        mo.image(os.path.join(public_root,f"data_exploration_{plate}_{exp}_gene_normalized_all_heatmap.png"))
    ]
    )
    return


@app.cell
def _(exp, mo, os, plate, public_root):
    mo.vstack([
        mo.md("### Heatmap - per guide"),
        mo.image(os.path.join(public_root,f"data_exploration_{plate}_{exp}_guide_normalized_all_heatmap.png"))
    ]
    )
    return


@app.cell
def _(exp, mo, os, pd, plate, public_root):
    gene_df = pd.read_csv(os.path.join(public_root,f"profiles_{plate}_{exp}_gene_normalized_feature_selected_byall.csv.gz"))
    mo.vstack([
        mo.md("### Per-Gene Profile Dataframe - explore!"),
        gene_df
    ]
    )
    return


@app.cell
def _(exp, mo, os, pd, plate, public_root):
    guide_df = pd.read_csv(os.path.join(public_root,f"profiles_{plate}_{exp}_guide_normalized_feature_selected_byall.csv.gz"))
    mo.vstack([
        mo.md("### Per-Guide Profile Dataframe - explore!"),
        guide_df
    ]
    )
    return


@app.cell
def _(exp, glob, os, plate):
    gene_radar_plot_list = glob.glob(f"**/data_exploration/gene_radar_plots/{exp}/{plate}/*_radar.png",recursive=True)
    gene_radar_plot_dict = {os.path.split(x)[1]:x for x in gene_radar_plot_list}
    gene_radar_plot_selectors = list(gene_radar_plot_dict.keys())
    gene_radar_plot_selectors.sort()
    return gene_radar_plot_dict, gene_radar_plot_selectors


@app.cell
def _(gene_radar_plot_selectors, mo):
    gene_radar_plot_el = mo.ui.dropdown(options=gene_radar_plot_selectors, value = gene_radar_plot_selectors[0], label = 'Choose a gene-level radar plot to view')
    mo.vstack([
        mo.md("### Per-Gene Radar Plots - pick a gene and see what has changed!"),
        gene_radar_plot_el
    ]
    )
    return (gene_radar_plot_el,)


@app.cell
def _(gene_radar_plot_dict, gene_radar_plot_el, mo):
    mo.image(gene_radar_plot_dict[gene_radar_plot_el.value])
    return


@app.cell
def _(exp, glob, os, plate):
    guide_radar_plot_list = glob.glob(f"**/data_exploration/guide_radar_plots/{exp}/{plate}/*_radar.png",recursive=True)
    guide_radar_plot_dict = {os.path.split(x)[1]:x for x in guide_radar_plot_list}
    guide_radar_plot_selectors = list(guide_radar_plot_dict.keys())
    guide_radar_plot_selectors.sort()
    return guide_radar_plot_dict, guide_radar_plot_selectors


@app.cell
def _(guide_radar_plot_selectors, mo):
    guide_radar_plot_el = mo.ui.dropdown(options=guide_radar_plot_selectors, value = guide_radar_plot_selectors[0], label = 'Choose a guide-level radar plot to view')
    mo.vstack([
        mo.md("### Per-guide Radar Plots - pick a guide and see what has changed!"),
        guide_radar_plot_el
    ]
    )
    return (guide_radar_plot_el,)


@app.cell
def _(guide_radar_plot_dict, guide_radar_plot_el, mo):
    mo.image(guide_radar_plot_dict[guide_radar_plot_el.value])
    return


@app.cell
def _(mo):
    mo.md("""
    ## Selected measurements of barcode quality calling
    ### Raw stats for this sub-experiment
    """)
    return


@app.cell
def _(exp, json, os, public_root):
    data_stats = json.load(open(os.path.join(public_root,f"summary_data_Data_Stats_{exp}.json")))
    data_stats
    return


@app.cell
def _(mo, os, public_root):
    mo.vstack([
        mo.md("### Spot quality bar graph"),
        mo.image(os.path.join(public_root,f"summary_figures_Quality_Spots_by_data_split.png"))
    ]
    )
    return


@app.cell
def _(mo, os, public_root):
    mo.vstack([
        mo.md("### Cell quality bar graph"),
        mo.image(os.path.join(public_root,f"summary_figures_Quality_Cells_by_data_split.png"))
    ]
    )
    return


@app.cell
def _(mo, os, plate, public_root):
    mo.vstack([
        mo.md("### Bias in barcoding errors - across all experiments"),
        mo.image(os.path.join(public_root,f"13_preprocessing_{plate}_observed_minus_exp_nucleotide_frequency.png"))
    ]
    )
    return


@app.cell
def _(mo):
    mo.md("""
    ## Selected measurements of overall experiment quality
    """)
    return


@app.cell
def _(batch, mo, os, plate, public_root):
    mo.vstack([
        mo.md("### Phenotyping illum corr (not experiment specific)"),
        mo.image(os.path.join(public_root,f"1_illumination_painting_{plate}_painting.{batch}_{plate}.montage.png"))
    ]
    )
    return


@app.cell
def _(batch, mo, os, plate, public_root):
    mo.vstack([
        mo.md("### Barcoding illum corr (not experiment specific)"),
        mo.image(os.path.join(public_root,f"5_illumination_barcoding_{plate}_barcoding.{batch}_{plate}.montage.png"))
    ]
    )
    return


@app.cell
def _(mo, well_list):
    barcode_alignment_plot_el = mo.ui.dropdown(options=well_list, value = well_list[0], label = 'Choose a well whose alignment you want to view; also used for the cross-arm alignment below')
    mo.vstack([
        mo.md("### Barcoding alignment plots - Cycles 2-12 to 1"),
        barcode_alignment_plot_el
    ]
    )
    return (barcode_alignment_plot_el,)


@app.cell
def _(barcode_alignment_plot_el, glob, mo, plate):
    barcode_alignment_plot_list = glob.glob(f"**/8_alignment_barcoding_{plate}*{barcode_alignment_plot_el.value}_*.png",recursive=True)
    barcode_alignment_plot_list.sort()
    mo.vstack([
    mo.hstack([mo.image(barcode_alignment_plot_list[0],width=250),mo.image(barcode_alignment_plot_list[1],width=250),mo.image(barcode_alignment_plot_list[2],width=250)]),
    mo.hstack([mo.image(barcode_alignment_plot_list[3],width=250),mo.image(barcode_alignment_plot_list[4],width=250),mo.image(barcode_alignment_plot_list[5],width=250)]),
    mo.hstack([mo.image(barcode_alignment_plot_list[6],width=250),mo.image(barcode_alignment_plot_list[7],width=250),mo.image(barcode_alignment_plot_list[8],width=250)]),
    mo.hstack([mo.image(barcode_alignment_plot_list[9],width=250),mo.image(barcode_alignment_plot_list[10],width=250)]),
            ]
        )
    #mo.image(guide_radar_plot_dict[guide_radar_plot_el.value])
    return


@app.cell
def _(barcode_alignment_plot_el, mo, os, plate, public_root):

    mo.vstack([
        mo.md("### Cross-arm alignment plots"),
    mo.image(os.path.join(public_root,f"9_alignment_combined_{plate}_{plate}-{barcode_alignment_plot_el.value}_AlignCombined_QC.png"))
    ]
    )
    return


@app.cell
def _(exp, mo, os, public_root):
    mo.vstack([
        mo.md("### Cell count per site and well"),
    mo.image(os.path.join(public_root,f"summary_figures_plate_layout_cells_count_{exp}.png"),width=150)
    ]
    )
    return


@app.cell
def _(exp, mo, os, public_root):
    mo.vstack([
        mo.md("### Over-confluency detection per site and well"),
    mo.image(os.path.join(public_root,f"summary_figures_plate_layout_percent_confluent_{exp}.png"),width=150)
    ]
    )
    return


@app.cell
def _(exp, mo, os, public_root):
    mo.vstack([
        mo.md("### Percent of cells with no spots per site and well"),
    mo.image(os.path.join(public_root,f"summary_figures_plate_layout_PercentEmpty_{exp}.png"),width=150)
    ]
    )
    return


@app.cell
def _(glob, mo):
    image_list = glob.glob("scripts/**/*.png",recursive=True)
    generic_image_picker = mo.ui.dropdown(options=image_list, value = image_list[0], label = "Pick an image to view")
    mo.vstack([
        mo.md(f"## Want to exlore more? Choose any image generated across all experiments ({len(image_list)} in total)"),
        generic_image_picker
    ]
    )
    return (generic_image_picker,)


@app.cell
def _(generic_image_picker, mo):
    mo.image(generic_image_picker.value)
    return


if __name__ == "__main__":
    app.run()
