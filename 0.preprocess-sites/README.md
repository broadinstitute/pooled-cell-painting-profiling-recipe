# Preprocessing Sites

In this first module of the recipe, we perform a series of steps to prepare data for profile generation.
The `preprocessing_config.yaml` controls all the pertinent details for this module.

<p align="center">
<img src="https://raw.githubusercontent.com/broadinstitute/pooled-cell-painting-profiling-recipe/82eaf532e7a3ab145c4b821268c13c531b693dcb/0.preprocess-sites/media/preprocessing_workflow.png" width="500">
</p>

## Step 0 - Prefilter Features

Many CellProfiler features per cell are not true morphology features.
This step flags all of these features for downstream removal.

## Step 1 - Process Spots

This step will process all in situ sequencing (ISS) spots, corresponding to CRISPR barcodes, across all imaging sites.
We output several figures describing barcode quality, barcode counts, cell quality, and cell assignment counts.
As well, we also output several summary text files to be visualized in a later step.

## Step 2 - Process Cells

This step will process all CellProfiler single cell features across all imaging sites.
We output two text files including cell metadata and cell count summary.

## Step 3 - Visualize Summary

This step compiles all the results from the spot and cell processing and outputs a series of visualizations providing an experiment-wide view of cell counts and cell quality.

## Step 4 - Image and Segmentation Quality Control

This step outputs a series of QC metrics pertaining to the actual image acquisition and image segmentation quality.
The QC metrics are visualized by site across wells per experiment.
