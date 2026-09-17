# Pooled Morphological Profiling Workflow

A step-by-step data processing pipeline for pooled morphological profiling (including Cell Painting) data.

## Overview

### Inputs

This worfkflow inputs morphological measurements from single cell objects and in situ sequencing (ISS, also referred to as sequencing by synthesis or SBS) measurements/barcode calls from a pooled screens.
Measurements are made by CellProfiler software using a [template pipeline](https://github.com/broadinstitute/pooled-cell-painting-image-processing/tree/master/pipelines/12cycles) (that can be modified as needed) including custom CellProfiler plugin modules for [channel balancing](https://github.com/CellProfiler/CellProfiler-plugins/blob/master/active_plugins/compensatecolors.py) and [barcode calling](https://github.com/CellProfiler/CellProfiler-plugins/blob/master/active_plugins/callbarcodes.py).

Example datasets for use with this workflow can be found on the [Cell Painting Gallery](https://broadinstitute.github.io/cellpainting-gallery/overview.html) and include `cpg0021-periscope` and `cpg0032-pooled-rare`.

### Recipe Steps

The workflow consists of six steps:

0. [Image QC](scripts/0.image-qc.py).
Image metrics from CellProfiler outputs are used to generate QC plots and reports.
Includes checking for cell confluency, image focus, image saturation, and phenotyping images to genotyping images alignment.
1. [Processing of SBS data](scripts/1.process-SBS.py).
Reads in SBS foci data from CellProfiler outputs and generates per-cell barcode/gene assignment and per-site SBS metrics.
2. [Merging SBS data and morphology data](scripts/2.merge-single-cells.py).
Reads in phenotyping data from CellProfiler outputs and merges it with SBS data to generate a single-cell level dataset.
Outputs a list of folders/sites that were used for single-cell dataset so that subsequent steps do not need access to the CellProfiler output data.
3. [Summarizing SBS data](scripts/3.summarize-SBS.py).
Creates QC plots and reports for SBS data and cell assignment.
Generates a total barcode count summary for comparison to NGS.
4. [Aggregating single cell profiles](scripts/4.aggregate.py).
Creates guide- and gene-level profiles by aggregating single-cell profiles using Pycytominer.
5. [Normalizing profiles](scripts/5.normalize.py).
Normalizes all profiles on a per-plate basis using Pycytominer.
6. [Feature selection](scripts/6.feature-selection.py).
Performs feature selection on normalized profiles using Pycytominer.
7. [Explore](scripts/7.explore.py).
Performs preliminary exploration of biological relationships in the profiles output by the recipe.

### Logging

The recipe includes creation of a log file for each step.
The file is named after the step (e.g. 0.prefilter-features.log) and saves in a logs/ folder.
The file logs progress information, warnings, and uncaught exceptions.

## Running the Recipe

The recipe supports running the workflow using Python in a conda environment or using Nextflow in either a conda environment or in a Docker.

### Step 1: Initialize the computational environment

This step can be skipped if you are using Nextflow with Docker.

Install [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).
We use conda as an environment manager.

```bash
# Install computational environment
conda env create --force --file environment.yml

# Initialize the environment
conda activate pooled-profiling
```

### Step 2: Configure the workflow

Edit the [Experiment configuration file](config/experiment.json) to set parameters specific to your experiment.
More information about the parameters can be found in the [config documentation](config/docs/experiment_README.md).

Note that there is also a [Defaults configuration file](config/defaults.json) that contains parameters that are less likely to need changing between experiments.
These include default assumptions about the structure of this repository and CellProfiler naming conventions.

### Step 3: Run the workflow

#### Run the workflow with Python

`python run.py`

If you would like to pass configuration files that have custom location or naming, you can do so with the following command:

```python
python run.py --defaults-config path/to/defaults.json --experiment-config path/to/experiment.json
```

#### Run the workflow with Nextflow

##### Nextflow and Python

The same recipe steps can be run as a [Nextflow](https://www.nextflow.io/) pipeline (`main.nf`) instead of `run.py`.
Each step is still one of the `scripts/N.*.py` files above; Nextflow just orchestrates which steps run and in what order, based on the `perform_*` flags in your experiment config.

Install Nextflow (see the [Nextflow docs](https://www.nextflow.io/docs/latest/install.html)), activate the `pooled-profiling` conda environment from Step 1, then run:

```bash
nextflow run main.nf
```

To use custom configuration file locations:

```bash
nextflow run main.nf --defaults_config path/to/defaults.yaml --experiment_config path/to/experiment.json
```

##### Nextflow and Docker

Instead of installing the `pooled-profiling` conda environment locally (Step 1), you can run every step inside a Docker container built from the same `environment.yml`. Install [Docker](https://docs.docker.com/get-docker/) and make sure the daemon is running, then build the image once:

```bash
docker build -t pooled-profiling .
```

Then run the workflow with the `docker` profile instead of activating a conda environment:

```bash
nextflow run main.nf -profile docker
```

The container only contains the Python environment; your repo checkout and data directories (`file_location`, `library_location`) are bind-mounted in automatically based on your config files, so editing scripts or configs never requires rebuilding the image — only changes to `environment.yml` do.
