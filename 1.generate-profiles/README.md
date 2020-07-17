# Generate Profiles

In this second module of the recipe, we perform a series of steps to convert single cell data into morphology profiles.
The `profiling_config.yaml` controls all the pertinent details for this module.

<p align="center">
<img src="https://raw.githubusercontent.com/broadinstitute/pooled-cell-painting-profiling-recipe/82eaf532e7a3ab145c4b821268c13c531b693dcb/1.generate-profiles/media/profiling_workflow.png" width="500">
</p>

## Step 0 - Merge Single Cells

This step merges compartment data and metadata, filters cells and features, and outputs single cell morphology profiles.

## Step 1 - Aggregate

Combine single cell data by averaging morphology features to form "aggregate" profiles.
We typically form aggregate profiles for both gene- and guide-level perturbations.

## Step 2 - Normalize

Perform a normalization scheme for the provided data level.
The data levels include single cell, genes, and guides.

## Step 3 - Feature Select

This step performs a series of feature selection steps to isolate the most pertinent morphology features in the data.
