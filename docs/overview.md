# Overview

A step-by-step data processing pipeline for pooled morphological profiling (including Pooled Cell Painting) data.

## Inputs

The workflow takes morphological measurements of single cell objects and in situ sequencing (ISS, also referred to as sequencing by synthesis or SBS) measurements/barcode calls from pooled screens.
Measurements are made by CellProfiler using a [template pipeline](https://github.com/broadinstitute/pooled-cell-painting-image-processing/tree/master/pipelines/12cycles) (that can be modified as needed) including custom CellProfiler plugin modules for [channel balancing](https://github.com/CellProfiler/CellProfiler-plugins/blob/master/active_plugins/compensatecolors.py) (for 4-color chemistry) and [barcode calling](https://github.com/CellProfiler/CellProfiler-plugins/blob/master/active_plugins/callbarcodes.py) (for 2-, 3-, or 4- color chemistry).

The CellProfiler output is expected to be nested as `file_location/BATCH/PLATE-WELL-SITE/*.csv`, with an `Image.csv`, `Foci.csv` (and, if barcode matching was done during image analysis, `BarcodeFoci.csv`), and one .csv per compartment (e.g. `Cells.csv`, `Nuclei.csv`, `Cytoplasm.csv`) in each site folder.

Example datasets for use with this workflow can be found on the [Cell Painting Gallery](https://broadinstitute.github.io/cellpainting-gallery/overview.html) and include `cpg0021-periscope` and `cpg0032-pooled-rare`.

## Recipe steps

| Step | Script | Purpose |
| --- | --- | --- |
| 0 | [Image QC](0.image-qc.md) | QC plots and reports from CellProfiler image metrics |
| 1 | [Process SBS](1.process-SBS.md) | Per-cell barcode/gene assignment and per-site SBS metrics |
| 2 | [Merge single cells](2.merge-single-cells.md) | Merge SBS and morphology data into single-cell datasets |
| 3 | [Summarize SBS](3.summarize-SBS.md) | QC plots and reports for SBS data and cell assignment |
| 4 | [Aggregate](4.aggregate.md) | Guide- and gene-level profiles |
| 5 | [Normalize](5.normalize.md) | Per-plate normalization |
| 6 | [Feature select](6.feature-select.md) | Remove uninformative features |
| 7 | [Explore](7.explore.md) | UMAPs, heatmaps, and mAP of the final profiles |

Each step is configured with the [configuration files](configs.md) and can be turned on or off with the `perform_*` flags in `experiment.json`.
Each step reads outputs of the preceding step(s), so steps can be re-run individually.
Step 2 is the last step that reads from the CellProfiler output folders.

## Output folder structure

Outputs are written under `directory_structure: root` (default `data/`) in `defaults.yaml`:

```text
data/
├── SBS/                 step 1 (and 2): per-site processed SBS data
├── single_cell/         steps 0, 2, 5: single-cell data, alignments
├── summary_data/        steps 0-3: tables and lists of sites
├── summary_figures/     steps 0, 3: QC figures
├── profiles/            steps 4-6: aggregated, normalized, feature selected profiles
└── data_exploration/    step 7: UMAPs, heatmaps, mAP
```

## Logging

The recipe creates a log file for each step.
The file is named after the step (e.g. `1.process-SBS.log`) and saved in a `logs/` folder.
The file logs progress information, warnings, and uncaught exceptions.

## Allowed skips

Steps 0-4 tolerate sites that fail to load.
Each of these steps has an `allowed_skips` setting in `defaults.yaml`; once that many sites have been skipped because of errors the step stops with an error.
