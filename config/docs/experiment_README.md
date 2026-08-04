# Documentation: experiment.json Configuration

Detailed information on how to customize the `experiment.json` for each Pooled morphological profiling experiment.

`file_location`: This workflow assumes that input data for the workflow was output from image analysis with a folder per site in this parent folder.
e.g. `"/Users/eweisbar/Desktop/demo/analysis"`

`control_genes`: A list of strings corresponding to the gene name given to controls in your experiment.
These are used in the normalization step to determine which cells to normalize against.
e.g. `["NT","nontargeting"]`

`drop_barcodes`: A list of any barcodes that should be excluded from analysis.
They should be the length of the barcodes called in the experiment.
e.g. `["AAAATTTTCCCCGGGG"]`

`data_sets`: A dictionary defining the different data subsets to be processed in this experiment.
If all of your data should be processed together, you should define a single data set that includes all batches, plates, and wells.
If you would like to process different parts of your data separately (for example, different wells received a different drug co-treatment), you can define multiple data sets here. \
e.g. single set: `{"Experiment1":{"Batch1":{"Plate1":["Well1","Well2","Well3","Well4","Well5","Well6"]}}}` \
e.g. multiple sets: `{"Untreated":{"20210422_6W_CP257":{"CP257A":["Well1","Well2"]}},"Drug_Treated":{"20210422_6W_CP257":{"CP257A":["Well3","Well4"]}}}`

`compartments`: A list of the cellular compartments that were measured in CellProfiler and should be included in the profiling.
The first compartment in the list must be the experimental unit (usually "Cells") and must be the parent object for the SBS Foci.
All compartments must have a Parent or Child relationship to the experimental unit so that cell_id_cols (default ImageNumber and ObjectNumber) can be used to merge them.
If the Parent/Child relationship is not inherent in the CellProfiler output through IdentifyPrimary/Secondary/Tertiary object modules, you will need to be sure you included a RelateObjects module to relate them back to the experimental unit. \
e.g. standard Cell Painting `["Cells","Nuclei","Cytoplasm"]` \
e.g. with additional organelle `["Cells","Nuclei","Cytoplasm","Mitochondria"]`

`cell_quality_method`: The method to use for determing the quality of cell assigment.
Cell quality methods are defined in [utils/cell_quality_utils.py](../utils/cell_quality_utils.py).
e.g. `"simple_T7"`

`overwrite_files`: Should existing output files be overwritten?
Note that this is only respected by `1.process-SBS.py` and `2.merge-single-cells.py`.
e.g. `true` or `false`

`infer_empty_sites`: Should the workflow attempt to infer empty sites in the output?
Workflows that perform square crops of a circular acquistion and/or workflows where the morphology acquisition map is different from the SBS acquisition map often have sites that completely black.
When run through the standard CellProfiler analysis pipeline, these will produce Image.csv's but not compartment .csv's.
If `true`, the pipeline will not count these sites as errors in the `allowed skips` limit.
e.g. `true` or `false`

`perform_step`: Should the step be performed?
This allows skipping of steps that have already been completed.
e.g. `true` or `false`
