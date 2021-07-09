import pathlib
import pandas as pd

from io_utils import read_csvs_with_chunksize

def load_single_cell_compartment_csv(compartment_dir, compartment, metadata_cols):
    """
    Load and process columns for CellProfiler output data

    Arguments:
    compartment_dir - path location of where the compartment csv files are stored
    compartment - string representing the compartment to load (e.g. cytoplasm)
    metadata_cols - a list of columns to add `Metadata_` prefix.
        Note the entries should not already be prefixed by compartment
        (e.g. AreaShape and not Cells_AreaShape)

    Output:
    A compartment dataframe with compartment prefixed column names
    """
    # Setup compartment file
    compartment = compartment.capitalize()
    compartment_file = pathlib.Path(compartment_dir, f"{compartment}.csv")

    # Load compartment data
    compartment_df = read_csvs_with_chunksize(compartment_file)
    compartment_df.columns = [f"{compartment}_{x}" for x in compartment_df.columns]

    # Identify and rename metadata_cols
    metadata_rename = {}
    for col in metadata_cols:
        metadata_col = f"Metadata_{compartment}_{col}"
        metadata_rename[f"{compartment}_{col}"] = metadata_col

    compartment_df = compartment_df.rename(metadata_rename, axis="columns")

    return compartment_df


def merge_single_cell_compartments(compartment_df_dict, merge_info_dict, id_cols):
    """
    Merge single cell compartment csvs using specified merge columns

    Arguments:
    compartment_df_dict - dictionary of pandas dataframes keyed by compartment name
    merge_info_dict - stores merge information and is loaded directly from config file
    id_cols - list of columns to identify single cell compartments, loaded from config

    Output:
    A single merged dataframe of all single cell measurements across compartments
    """
    # Extract expected merge information
    link_compartment = merge_info_dict["linking_compartment"].capitalize()
    linking_columns = merge_info_dict["linking_columns"]
    image_col = merge_info_dict["image_column"]
    link_compartment_col = f"Metadata_{link_compartment}_{image_col}"
    linker_df = compartment_df_dict[link_compartment]

    # Perform the merge given the linking compartment columns
    for compartment, compartment_link in linking_columns.items():
        # Pull compartment dataframe to link from the given dictionary
        to_link_df = compartment_df_dict[compartment.capitalize()]

        # Setup the columns to match between dataframes
        link_merge_cols = (link_compartment_col, compartment_link)
        to_link_compartment_cols = [
            f"Metadata_{compartment.capitalize()}_{x}" for x in id_cols
        ]

        # Merge compartment dataframes
        linker_df = linker_df.merge(
            to_link_df,
            left_on=link_merge_cols,
            right_on=to_link_compartment_cols,
            how="inner",
        )

    return linker_df
