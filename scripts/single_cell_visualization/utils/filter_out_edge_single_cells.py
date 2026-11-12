"""
@author: mhaghigh

* edge cells:
  - Are the cells which their nuclei center is less than a specified margin,
    - The margin could be an optional input or if not specified, 
    - By default, it is the 90th percentile of the "Cells_AreaShape_MajorAxisLength" feature devided by 2

"""

import pandas as pd
import numpy as np

######################## function to remove cells on the border
def edgeCellFilter(df_input, location_column_X, location_column_Y, image_width=None, edge_margin=None):
    """
    remove cells close to the image borders
    edge cells are the cells which their nuclei center is less than a specified margin,

    Inputs:
        df_input
        image_width (optional) if not given as the input it would be inferred from the input dataframe columns
        edge_margin (optional) if not given as the input, it would be calculated
                from the input dataframe columns
                - By default, it is the 90th percentile of the "Cells_AreaShape_MajorAxisLength" feature devided by 2

    """

    if not image_width:
        metadata_cols_4size = df_input.columns[
            df_input.columns.str.contains("Width|ImageSize")
        ]
        if len(metadata_cols_4size):
            image_width = df_input[metadata_cols_4size[0]].values[0]
        else:
            raise Exception(
                "No metadata columns for inferring image size are detected! Please enter edge_margin as an input!"
            )

    if not edge_margin:
        if "Cells_AreaShape_MajorAxisLength" in df_input.columns:
            edge_margin = int(
                np.percentile(df_input.Cells_AreaShape_MajorAxisLength.values, 80) / 2
            )
            print (f"Dropping cells within {edge_margin} pixels from the image edges")
        else:
            raise Exception(
                "The Cells_AreaShape_MajorAxisLength column is not detected! Please enter image_width as an input!"
            )

    df_filtered_edge_cells = df_input.loc[
        ~(
            (df_input[location_column_X] > (image_width - edge_margin))
            | (df_input[location_column_X] < (edge_margin))
            | (df_input[location_column_Y] > (image_width - edge_margin))
            | (df_input[location_column_Y] < (edge_margin))
        ),
        :,
    ].reset_index(drop=True)
    print (f"Dropped {len(df_input)-len(df_filtered_edge_cells)} edge cells out of {len(df_input)} total cells that were within {edge_margin} pixels from the image edges")

    return df_filtered_edge_cells, edge_margin
