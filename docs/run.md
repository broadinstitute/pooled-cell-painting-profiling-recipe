# Running the recipe

The recipe supports running the workflow using Python in a conda environment or using Nextflow in either a conda environment or in a Docker container.

## Step 1: Initialize the computational environment

This step can be skipped if you are using Nextflow with Docker.

Install [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).
We use conda as an environment manager.

```bash
# Install computational environment
conda env create --force --file environment.yml

# Initialize the environment
conda activate pooled-profiling
```

## Step 2: Configure the workflow

Edit `config/experiment.json` to set parameters specific to your experiment.
`config/defaults.yaml` contains parameters that are less likely to need changing between experiments.
See [Configuration](configs.md) for every parameter.

## Step 3: Run the workflow

### Run the workflow with Python

```bash
python run.py
```

To pass configuration files that have custom location or naming:

```bash
python run.py --defaults-config path/to/defaults.yaml --experiment-config path/to/experiment.json
```

`run.py` checks that `experiment.json` contains all required keys and then runs each step in order for which the matching `perform_*` flag is `true`.
It stops if any step exits with an error.

A single step can also be run directly by passing it the two config paths:

```bash
python scripts/1.process-SBS.py config/defaults.yaml config/experiment.json
```

### Run the workflow with Nextflow

The same steps can be run as a [Nextflow](https://www.nextflow.io/) pipeline (`main.nf`) instead of `run.py`.
Each step is still one of the `scripts/N.*.py` files; Nextflow orchestrates which steps run and in what order, based on the `perform_*` flags in your experiment config.

Install Nextflow (see the [Nextflow docs](https://www.nextflow.io/docs/latest/install.html)), activate the `pooled-profiling` conda environment from Step 1, then run:

```bash
nextflow run main.nf
```

To use custom configuration file locations:

```bash
nextflow run main.nf --defaults_config path/to/defaults.yaml --experiment_config path/to/experiment.json
```

### Run Nextflow with Docker

Instead of installing the conda environment locally, you can run every step inside a Docker container built from the same `environment.yml`.
Install [Docker](https://docs.docker.com/get-docker/) and make sure the daemon is running, then build the image once:

```bash
docker build -t pooled-profiling .
```

Then run the workflow with the `docker` profile:

```bash
nextflow run main.nf -profile docker
```

The container only contains the Python environment; your repo checkout and data directories (`file_location`, `library_location`) are bind-mounted automatically based on your config files, so editing scripts or configs never requires rebuilding the image — only changes to `environment.yml` do.
