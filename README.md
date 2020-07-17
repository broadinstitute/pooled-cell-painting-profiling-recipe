# Pooled Cell Painting - Image-based Profiling Pipeline Recipe :woman_cook: :man_cook:

A step-by-step data processing pipeline for Pooled Cell Painting data.

## Ingredients

Data are the primary ingredients of science.
Here, our data come from a Pooled Cell Painting experiment.

In these experiments, the data are thousands of `.csv` files storing metadata and morphology measurements from millions of single cells.

There are two fundamental kinds of data ingredients:

1. Cells
2. Spots

The `Cells` ingredients represent morphology measurements for various cellular compartments for each segmented single cell.
The `Spots` ingredients represent in situ sequencing (ISS) results used for "cell calling".
Cell calling is the procedure that assigns a specific CRISPR barcode to an individual single cell.
Because the experiment is "pooled", there are thousands of CRISPR barcodes present in a single well.

These measurements for both data ingredients are currently made by CellProfiler software (using customized Pooled Cell Painting plugins).

## Recipe Steps

All cookbooks also include specific instructions, or steps, for each recipe.

Our recipe includes two modules:

1. [Preprocessing](0.preprocess-sites/)
2. [Profile generation](1.generate-profiles/)

The output data are structured in a way that includes measurements from many individual "sites" across a single plate.
Each site can be thought of as a single field of view that consists of many different images from the five Cell Painting channels, and four ISS channels across `n` cycles.
The number of cycles is determined as part of experimental design and is typically selected to ensure zero collisions between CRISPR barcodes.

The recipe steps first preprocess spots and cells, output quality control (QC) metrics, and perform filtering.
Next, in profile generation, single cell profiles are merged, aggregated, normalized and feature selected.
The final output of the pipeline are QC metrics, summary figures, and morphology profiles for each CRISPR guide.
These profiles will be used in downstream analyses for biological discovery.

## Usage

This recipe is designed to be used as a critical component of a Data Pipeline Welding procedure.

More specifically, this recipe will be linked together, via a [GitHub submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules), to a Pooled Cell Painting data repository.
The data repositories will be derived from the [Pooled Cell Painting template](https://github.com/broadinstitute/pooled-cell-painting-profiling-template).

More usage instructions can be found in the template repo linked above.
Briefly, the goal of the weld is to tightly couple the Pooled Cell Painting processed data to versioned code that performed the processing.
This recipe is the versioned code and a GitHub submodule links the recipe by commit hash.

The recipe is interacted with via a series of configuration yaml files defined in the data repository.
