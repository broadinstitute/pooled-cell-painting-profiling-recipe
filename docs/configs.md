# Configuration

The recipe is configured with two files:

- `config/experiment.json` holds parameters that change with every experiment (where the data is, which wells to process, which genes are controls, which steps to run).
- `config/defaults.yaml` holds parameters that are less likely to change between experiments, such as output folder names, CellProfiler column names, and per-step options.

Pass alternative files with `--defaults-config` and `--experiment-config` (see [Running the recipe](run.md)).
`config/experiment.json` keys beginning with `_commenton_` are in-file comments and are ignored.

## experiment.json

### Input data

`file_location`: The parent folder for a set of .csv's that have been output into nested `/BATCH/PLATE-WELL-SITE` within the parent folder. \
e.g. `"/Users/me/Desktop/demo/analysis"`

`external_foci_location`: If you have made measurements on foci that did not happen during the same image_analysis step but instead are brought in from another step, this is the parent folder for those .csv's.
Note that the Foci objects must match 1:1 by Object Number to use `external_foci_location`.
Supports both the same `/BATCH/PLATE-WELL-SITE` nesting as `file_location`, or `/PLATE/PLATE-WELL-SITE` nesting.
Files may have a prefix before `Foci.csv` (e.g. `BarcodePreprocessing_Foci.csv`). `*BarcodeFoci.csv` is ignored when locating the foci file.
Leave blank if your foci measurements are in the same parent folders as other measurements. \
e.g. `"/Users/me/Desktop/separate_foci_folder/"`

`data_sets`: A dictionary defining the different data subsets to be processed in this experiment, in the form `{SETNAME: {BATCH: {PLATE: [WELL1, WELL2]}}}`.
Each data set must contain exactly one batch but can contain any number of plates and wells.
If all of your data should be processed together, define a single data set that includes all plates and wells.
If you would like to process different parts of your data separately (for example, different wells received a different drug co-treatment), define multiple data sets. \
e.g. single set: `{"Experiment1":{"Batch1":{"Plate1":["A01","A02","A03"]}}}` \
e.g. multiple sets: `{"Untreated":{"20210422_6W_CP257":{"CP257A":["Well1","Well2"]}},"Drug_Treated":{"20210422_6W_CP257":{"CP257A":["Well3","Well4"]}}}`

`compartments`: A list of the cellular compartments that were measured in CellProfiler and should be included in the profiling.
The first compartment in the list must be the experimental unit (usually `Cells`) and must be the parent object for the SBS Foci.
All compartments must have a Parent or Child relationship to the experimental unit so that `cell_id_cols` (default `ImageNumber` and `ObjectNumber`) can be used to merge them.
If the Parent/Child relationship is not inherent in the CellProfiler output through IdentifyPrimary/Secondary/Tertiary object modules, you will need to be sure you included a RelateObjects module in the image analysis pipeline to relate them back to the experimental unit. \
e.g. standard Cell Painting `["Cells","Nuclei","Cytoplasm"]` \
e.g. with additional organelle `["Cells","Nuclei","Cytoplasm","Mitochondria"]`

### Controls

`control_genes`: Required. A list of strings corresponding to the gene name given to negative controls in your experiment.
These are used in the normalization step to determine which cells to normalize against (when `by_samples` is `control_genes`) and are highlighted in the exploration step. \
e.g. `["NT","nontargeting"]`

`poscon_genes`: Optional (may be an empty list `[]`, but the key must be present). A list of gene names to treat as positive controls.
Used in the exploration step to highlight positive controls and calculate similarity to them. \
e.g. `["POLR2A"]`

`drop_barcodes`: A list of any barcodes that should be excluded from analysis.
They should be the length of the barcodes called in the experiment.
Not supported when `library_structure` has more than one segment (dropping is skipped). \
e.g. `["AAAATTTTCCCCGGGG"]`

### Barcode calling and cell quality

`match_to_library`: `true` if barcode calls need to be matched to a barcode library in the recipe (during 1.process-SBS); `false` if matching already happened during image analysis (a `BarcodeFoci.csv` is then read for each site).

`library_structure`: A dictionary where keys are library columns to match, in order, and values are `{"n": length to match, "namecol": column with names to match}`.
Barcode calls are split into consecutive segments of length `n`, one per key.
Use a single key for a single-barcode library, or multiple keys for libraries with multiple barcodes per construct (e.g. dual guides). \
e.g. `{"iBAR_2":{"n":7,"namecol":"target_symbol"},"iBAR_1":{"n":5,"namecol":"target_symbol"}}` \
e.g. `{"Barcode":{"n":12,"namecol":"Gene"}}}`

`library_location`: Path to the barcode library .csv. Must contain every key of `library_structure` as a column, as well as each `namecol`.

`cell_quality_method`: The method used to determine the quality of each cell's barcode assignment.
Options are `simple`, `simple_T7`, `simple_plus`, and `multiples_T7`.
`multiples_T7` is the only method compatible with multiple-segment libraries.
Methods and the categories they produce are defined in [utils/cell_quality_utils.py](https://github.com/broadinstitute/pooled-cell-painting-profiling-recipe/blob/master/utils/cell_quality_utils.py). \
e.g. `"multiples_T7"`

`keep_cell_qualities`: A list of the cell quality category names (as produced by `cell_quality_method`) to keep in the single-cell data in step 2.
Cells in other categories are dropped before aggregation. \
e.g. for `multiples_T7`: `["Perfect","Good"]`; for `simple_T7`: `["Perfect","Great"]`

### Running

`overwrite_files`: Should existing output files be overwritten?
Note that this is only respected by `1.process-SBS.py` and `2.merge-single-cells.py`. \
e.g. `true` or `false`

`perform_process_qc`, `perform_process_SBS`, `perform_process_single_cells`, `perform_summarize_SBS`, `perform_aggregate`, `perform_normalize`, `perform_feature_select`, `perform_explore`: Should steps 0 through 7 (respectively) be performed?
This allows skipping of steps that have already been completed. \
e.g. `true` or `false`

## defaults.yaml

Note that structure is very important in `.yaml` files.
When editing the `.yaml` make sure to maintain dashes and indentations.
(When information is added on the same line, it makes a dictionary where the value is a string.
When there are multiple lines, the value is a list with string elements.)

### directory_structure

The naming for the output directories created by the workflow.
All are created within `root`.

`root`: Parent folder for all outputs.

`SBS`: Per-site processed SBS data (step 1).

`phenotyping`: Reserved for phenotyping data.

`data`: Summary data tables (.csv/.tsv/.json) produced by QC and summary steps.

`figures`: Summary figures produced by QC and summary steps.

`profiles`: Aggregated, normalized, and feature selected profiles.

`single_cell`: Merged single-cell data.

`exploration`: Output of the exploration step.

### core

`cell_quality_column`: The name of the column that contains the string describing the cell quality.

`cell_quality_index`: The name of the column that contains the integer describing the cell quality.

`cell_id_cols`: The list of columns that uniquely identify each cell in the single cell data.

`compression`: Compression to use when creating .csv files.
For more information and options, see Pandas DataFrame.to_csv documentation.

`float_format`: Decimal precision to use in writing output files.
For more information and options, see Python string formatting documentation.
Default `"%.5g"` keeps 5 significant digits.

`ignore_files`: Ignore any files with these names.
List the complete file name.
Default is `.DS_Store` because Macs make hidden files with that name.

### process

#### qc

`allowed_skips`: Number of sites you allow to be skipped because of an error during preprocessing before considering the step to have failed.

`infer_empty_sites`: Should the workflow attempt to infer empty sites in the output?
Workflows that perform square crops of a circular acquisition and/or workflows where the morphology acquisition map is different from the SBS acquisition map often have sites that are completely black.
When run through the standard CellProfiler analysis pipeline, these will produce Image.csv's but not compartment .csv's.
If `true`, the pipeline will not count these sites as errors in the `allowed_skips` limit, and later steps will skip them. \
e.g. `true` or `false`

`stack_alignment_chans`: A list of the two channel names that should be compared between phenotyping and barcoding images to determine alignment quality.
Channel names must appear in the `Correlation_Correlation_` column of `Image.csv`.

#### process_SBS

`allowed_skips`: Number of sites you allow to be skipped because of an error during preprocessing before considering the step to have failed.

`barcode_col`: The name of the column that contains the barcode assigned to a given SBS focus.

`call_col`: The name of the column that contains the barcode as called from the sequencing images (before matching to the library).

`gene_col`: The name of the column that contains the gene assigned to a given SBS focus.

`SBS_score_col`: The column that contains the match score for the barcode assigned to a given SBS focus.

`location_cols`: A list of the columns that define the location of each SBS focus (generally X and Y coordinates).

`foci_cols`: A list of the columns that contain the barcode that was called/assigned to a given SBS focus (listed first) and its corresponding gene or ID.

`spot_quality_method`: The method used to categorize the quality of individual spots when the library has multiple segments.
Currently only `simple` is supported.

#### single_cell

`allowed_skips`: Number of sites you allow to be skipped because of an error during preprocessing before considering the step to have failed.

`flag_cols`: A list of strings.
Phenotyping data drops all columns that contain any of these strings.
Note columns containing these strings may be used for other data processing, but will not be considered morphology features.

`save_single_file`: Boolean indicating whether to save all single cell data for a plate in one file.
May have both this and `save_per_site_file` set to true.
Cannot have both this and `save_per_site_file` set to false.

`save_per_site_file`: Boolean indicating whether to save single cell data in separate files per site.
May have both this and `save_single_file` set to true.
Cannot have both this and `save_single_file` set to false.

#### summarize_SBS

`allowed_skips`: Number of sites you allow to be skipped because of an error during preprocessing before considering the step to have failed.

#### aggregate

`allowed_skips`: Number of sites you allow to be skipped because of an error during preprocessing before considering the step to have failed.

`operation`: A string that is passed to Pycytominer's aggregate function to indicate how the data is aggregated.
Currently only supports `mean` or `median`.
See Pycytominer documentation for more information.

`features`: A list of strings that is passed to Pycytominer's aggregate function to indicate the list of features that should be aggregated.
If set to the string `infer`, Pycytominer will infer the features list instead of needing an explicit list passed.
See Pycytominer documentation for more information.
Note that Pycytominer assumes standard Cell Painting compartments of `Cell`, `Cytoplasm`, and `Nucleus`.

#### normalize

`output_bygene`: Boolean indicating whether to additionally write the normalized profiles split into one file per gene.

`method`: The normalization method passed to Pycytominer's normalize function (e.g. `standardize`, `robustize`, `mad_robustize`).
See Pycytominer documentation for more information.

`by_samples`: Which samples are used to calculate normalization statistics.
`control_genes` uses only samples whose gene is in `control_genes` in `experiment.json`.
`all` uses all samples.

`features`: A list of strings that is passed to Pycytominer's normalize function to indicate the list of features that should be normalized.
If set to the string `infer`, Pycytominer will infer the features list.

#### feature_select

`group`: The grouping at which you would like to perform feature selection.
All plates within this group will be used to calculate feature selection metrics together.
Supports `plate` (each plate separately), `group` (all plates within a data set), or `all` (all plates in all data sets).

`operations`: The list of feature selection operations in Pycytominer that you would like to use.
See Pycytominer documentation for more information.

`use_samples`: A list of samples to provide operation on.
If set to the string `all`, all samples are used.
See Pycytominer documentation for more information.

`features`: A list of strings that is passed to Pycytominer's feature_select function to indicate the list of features that should be selected from.
If set to the string `infer`, Pycytominer will infer the features list instead of needing an explicit list passed.
See Pycytominer documentation for more information.
Note that Pycytominer assumes standard Cell Painting compartments of `Cell`, `Cytoplasm`, and `Nucleus`.

`na_cutoff`: The proportion of missing values in a column to tolerate before removing.

`corr_threshold`: Float between (0, 1) to exclude features if they have a correlation above.
