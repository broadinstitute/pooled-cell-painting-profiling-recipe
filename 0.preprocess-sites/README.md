# Preprocessing Sites

In this first module of the recipe, we perform a series of steps to prepare data for profile generation.
The `preprocessing_config.yaml` controls all the pertinent details for this module.

<p align="center">
<img src="https://raw.githubusercontent.com/broadinstitute/pooled-cell-painting-profiling-recipe/068c7eae79f56a50732fdf173902dc596c95ce3f/0.preprocess-sites/media/preprocessing_workflow.png" width="500">
</p>

## Step 0 - Prefilter Features

Many CellProfiler features per cell are not true morphology features.
This step flags all of these features for downstream removal.
It removes all features that have the string in them.
(e.g. `Count` will remove `Count_Children_Foci`, `Count_Children_BarcodeFoci`, etc.)

## Step 1 - Process Spots

This step processes all in situ sequencing (ISS) spots, corresponding to CRISPR barcodes, across all imaging sites.
It outputs several figures describing barcode quality, barcode counts, cell quality, and cell assignment counts.
It also outputs several summary text files to be visualized in a later step.

## Step 2 - Process Cells

This step combines all CellProfiler single cell features measured in all compartments across all imaging sites.
It outputs two text files including a cell metadata summary and cell counts by quality category summary.

## Step 3 - Visualize Summary

This step compiles all the results from the spot and cell processing and outputs a series of visualizations providing an experiment-wide view of cell counts and cell quality.

## Step 4 - Image and Segmentation Quality Control

This step outputs a series of QC metrics pertaining to the actual image acquisition and image segmentation quality.
The QC metrics are visualized by site across wells per experiment.
It also outputs lists of sites that trigger QC flags for confluence, image saturation, and image alignment.
