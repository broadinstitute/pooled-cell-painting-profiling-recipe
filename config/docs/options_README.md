# Documentation: defaults.yaml Configuration

Detailed information on how to customize the `defaults.yaml` for each Pooled Cell Painting experiment.  

For more information on `.yaml`, read the [YAML page on Wikipedia](https://en.wikipedia.org/wiki/YAML).

Structure is very important in `.yaml` files.
When editing the `.yaml` make sure to maintain dashes and indentations.
(When information is added on the same line, it makes a dictionary where the value is a string.
When there are multiple lines, the value is a list with string elements.)

## directory structure

The naming for the output directories created by the workflow.

## core

`cell_quality_column`: The name of the column that contains the string describing the cell quality.

`cell_quality_index`: The name of the column that contains the integer describing the cell quality.

`cell_id_cols`: The list of columns that uniquely identify each cell in the single cell data.

`compression`: Compression to use when creating .csv files.
For more information and options, see Pandas DataFrame.to_csv documentation.

`float_format`: Decimal precision to use in writing output files.
For more information and options, see Python string formatting documentation.
Default `"%.5g"` keeps 5 decimal places.

`ignore_files:` Ignore any files with these names.
List the complete file name.
Default is `.DS_Store` because Macs make hidden files with that name.

## process

### qc

`allowed_skips`: Number of sites you allow to be skipped because of an error during preprocessing before considering the step to have failed.

`stack_alignment_chans`: A list of the two channel names that should be compared between phenotyping and barcodign images to determine alignment quality.

### process_SBS

`allowed_skips`: Number of sites you allow to be skipped because of an error during preprocessing before considering the step to have failed.

`barcode_col`: The name of the column that contains the barcode assigned to a given SBS focus.

`gene_col`: The name of the column that contains the gene assigned to a given SBS focus.

`SBS_score_cols`: The column that contains the match score for the barcode assigned to a given SBS focus.

`location_cols`: A list of the columns that define the location of each SBS focus (generally X and Y coordinates).

`foci_cols`: A list of the columns that contain the barcode that was called/asigned to a given SBS focus (listed first) and its corresponding gene.

### single_cell

`allowed_skips`: Number of sites you allow to be skipped because of an error during preprocessing before considering the step to have failed.

`flag_cols`: A list of strings.
Phenotyping data drops all columns that contain any of these strings.
Note columns containing these strings may be used for other data processing, but will not be considered morphology features.

`save_single_file`: Boolean indicating whether to save all single cell data in one file.
May have both this and `save_per_site_file` set to true.
Cannot have both this and `save_per_site_file` set to false.

`save_per_site_file`: Boolean indicating whether to save single cell data in separate files per site.
May have both this and `save_single_file` set to true.
Cannot have both this and `save_single_file` set to false.

### summarize_SBS

`allowed_skips`: Number of sites you allow to be skipped because of an error during preprocessing before considering the step to have failed.

### aggregate

`allowed_skips`: Number of sites you allow to be skipped because of an error during preprocessing before considering the step to have failed.

`operation`: A string that is passed to Pycytominer's aggregate function to indicate how the data is aggregated.
Currently only supports `mean` or `median`.
See Pycytominer documentation for more information.

`features`: A list of strings that is passed to Pycytominer's aggregate function to indicate the list of features that should be aggregated.
If set to the string `infer`, Pycytominer will infer the features list instead of needing an explicit list passed.
See Pycytominer documentation for more information.
Note that Pycytominer assumes standard Cell Painting compartments of `Cell`, `Cytoplasm`, and `Nucleus`.

### feature_select

`group`: The grouping at which you would like to perform feature selection.
All plates within this group will be used to calculate feature selection metrics together.
Supports `plate`, `group`, or `all`.

`operations`: The list of feature selection operations in Pycytominer that you would like to use.
See Pycytominer documentation for more information.

`use_samples`: A list of samples to provide operation on.
If set to the string `all`, all samples are used.
See Pycytominer documentation for more information.

`features`: A list of strings that is passed to Pycytominer's feature_select function to indicate the list of features that should be aggregated.
If set to the string `infer`, Pycytominer will infer the features list instead of needing an explicit list passed.
See Pycytominer documentation for more information.
Note that Pycytominer assumes standard Cell Painting compartments of `Cell`, `Cytoplasm`, and `Nucleus`.

`na_cutoff`: The proportion of missing values in a column to tolerate before removing.

`corr_threshold`: Float between (0, 1) to exclude features if they have a correlation above.
