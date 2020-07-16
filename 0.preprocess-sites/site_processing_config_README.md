Detailed information on how to customize site_processing_config.yaml for each experiment.  

When editing the .yaml make sure to maintain dashes and indentations. (When information is added on the same line, it makes a dictionary where the value is a string. When there are multiple lines the value is a list with string elements.)

## master_config  
*These settings are passed to a separate config file to give it information about this config file. Do not edit these.*  

## core:  
*This workflow assumes that input data was output from image analysis with a folder per site in the parent folder structure: project_dir/project/batch/*

`batch:` *The full batch tag as found in s3://imaging-platform/projects/PROJECT/  
e.g.* 20190919_6W_CP074A

`project_dir:` *This is the master folder that contains your input data. It is one level up from the project folder.  
e.g.* /Users/gway/work/projects/


`categorize_cell_quality:` *The cell categorization method you would like to use. The methods are described in cell_quality_utils.py.  
e.g.* simple  
*e.g.* simple_plus

`sites_per_image_grid_side:` *The number of images making one side of the square  grid of images being processed. If there are 100 input images then:  
e.g.* 10  

`compartments:` *The cellular compartments measured in the cell painting experiment. Default cell painting compartments are:  
e.g.*  
    - Cells  
    - Nuclei  
    - Cytoplasm  
*Additional possible compartments such as those used for SABER are:  
e.g.*    
    - Mito

`parent_cols:`  
*Enter every parent object/s for every compartment listed in compartments above (you may need to add additional key pairs for additional compartments beyond the default cells and cytoplasm e.g. for additional compartments beyond default imaged by SABER). Parent refers to relationship in order of identification. Generally, nuclei are identified first so they do not have parents. The parent columns can be found by searching the Compartment.csv files for columns following the format "Parent_Compartment."*  
`cells:` *Generally, nuclei are used to identify cells as a secondary object so the only parent object a cell will have is nuclei.  
e.g.*   
    - Parent_Nuclei  
`cytoplasm:` *Cytoplasm is generally a tertiary object identified by subtracting nuclei from cells so its parent objects are nuclei and cells.  
e.g.*  
    - Parent_Nuclei  
    - Parent_Cells  
`spots:` *Spots are generally labeled on a per-cell basis.  
e.g.*  
    - Parent_Cells  

`id_cols:`*All columns necessary for parsing objects. Objects are identified and numbered on a per-image basis. Therefore, generally both image number and object number are needed to parse objects.  
e.g.*    
    - ImageNumber  
    - ObjectNumber

`control_barcodes:`*The name/s used for non-targeting sgRNA's.    
e.g.*  
    - NT

`ignore_files:`*Ignore any files with these names. List the complete file name. Default is .DS_Store because Macs make hidden files with that name.  
e.g.*  
    - .DS_Store

## prefilter:
`perform:` *Do you want to perform 0.prefilter-features step? Set to* true *or* false.

`example_site:` *The name of a single example folder containing data from a single site that is found in project_dir/project/batch.  
e.g.* CP074A_A1-Site_1  
*e.g.* 151B2-B1-45

`output_basedir:` *Name the folder where you would like this step to output files. e.g.* data

`flag_cols:` *Searches every feature name and if it finds these it flags them. Flagged files are ignored during 0.prefilter-features.py though they may be used later on.  
e.g.*  
    - Barcode  
    - Location  
    - Count  
    - Resize *(Resized images are used for visual output and not measurements.)*

## process-spots:
`perform:` *Do you want to perform 1.process-spots step? Set to* true *or* false.

`output_basedir:` *Name the folder where you would like this step to output files. e.g.* data

`barcode_cols:` *The column that contains the barcodes (or columns that contain the barcodes) assigned to a given spot.  
e.g.*  
    - Barcode_MatchedTo_Barcode

`gene_cols:` *The column that contains the gene (or columns that contain the genes) assigned to a given spot.  
e.g.*  
    - Barcode_MatchedTo_GeneCode

`location_cols:` *The columns that define the location of each spot (generally X and Y coordinates).  
e.g.*  
    - Location_Center_X  
    - Location_Center_Y

`spot_score_cols:` *The column that contains the match score for the barcode (or columns that contain the match score for the barcodes) assigned to a spot.  
e.g.*  
    - Barcode_MatchedTo_Score

`foci_cols:`*The columns that contain the top match for barcode and corresponding gene.
e.g.*  
    - Barcode_BarcodeCalled  
    - Barcode_MatchedTo_ID

`cell_filter:`*Select the categories of cells that you wish to use for analysis. Categories are defined in cell_quality_utils.py. Make sure that the options you set here are available in the categorize_cells method that you set above in core config.  
Typical e.g.*  
    - Perfect  
    - Great  
*Other possible options e.g.*  
    - Imperfect *(only in simple)*  
    - Imperfect-High *(only in simple simple_plus)*  
    - Imperfect-Low *(only in simple_plus)*  
    - Bad

## process-cells:
`perform:` *Do you want to perform 2.process-cells step? Set to* true *or* false.

`prefilter-features:` *Do you want to use the prefiltering that you created in step 0? Set to* true *or* false.

`output_basedir:` *Name the folder where you would like this step to output files. e.g.* data

`sort_col:` *Column used for sorting data numerically. e.g.* Metadata_Cells_ObjectNumber  

#### merge_columns:
*CellProfiler outputs measurements for each compartment in a separate .csv. In this step we merge them together. merge_columns are the columns used to merge the .csvs.*  
`image_column:` *The name of the column that describes the image number. e.g.* ImageNumber  
`linking_compartment:` *The compartment used to link other compartments. e.g.* Cytoplasm  
`linking_columns:`*Column names that link your parent compartment to other compartments. You will need to set a column name for each compartment set in core: compartments. e.g.*  
`cells:` Metadata_Cytoplasm_Parent_Cells  
`nuclei:` Metadata_Cytoplasm_Parent_Nuclei  
`metadata_linking_columns:` *Column names that link ... e.g.*  
    - Metadata_Foci_site  
    - Metadata_Cells_ObjectNumber  

#### metadata_merge_columns:
*These are the columns that are used to link foci data to compartment data.*   
`foci_cols:` *The columns that describe the image number and the cell number in the foci data. e.g.*  
    - Metadata_Foci_ImageNumber  
    - Metadata_Foci_Parent_Cells  
`cell_cols:`*The columns that describe the image number and the cell number in the cell data. e.g.*  
    - Metadata_Cells_ImageNumber  
    - Metadata_Cells_ObjectNumber  
`cell_quality_col:` *The column that contains the numerical cell quality category. e.g.* Metadata_Foci_Cell_Quality

`foci_site_col:`*The column in the foci data that contains the site number e.g.* Metadata_Foci_site

## summarize-cells:
`perform:` *Do you want to perform 2.process-cells step? Set to* true *or* false.

`output_resultsdir:` *Name the folder where you would like this step to output files. e.g.* results

`output_figuresdir:` *Name the folder where you would like this step to output figures. e.g.* figures

## summarize-plate:
`perform:` *Do you want to perform 2.process-cells step? Set to* true *or* false.

`correlation_threshold:` *Set the minimum correlation for two images to have to pass the correlation filter. e.g.* .2

`painting_image_names:` *The names of all input images used for Cell Painting as they are named in CellProfiler. You can find image names as part of any measurement column in Image.csv (e.g. columns starting with Width_).  
e.g.*  
    - Pre_ConA  
    - Pre_DAPI_Painting  
    - Pre_Mito  
    - Pre_Phalloidin  
    - Pre_WGA  

`barcoding_cyles:` *The number of SBS cycles run in the experiment. e.g.* 12

`barcoding_prefix:` *The prefix at the start of all input images used for SBS as they are named in CellProfiler. You can find image names as part of any measurement column in Image.csv (e.g. columns starting with Width_) e.g.* CorrCycle
