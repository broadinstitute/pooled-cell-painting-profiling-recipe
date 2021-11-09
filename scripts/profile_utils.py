import pandas as pd


def sanitize_gene_col(metadata_df, gene_col, control_barcodes):
    """
    Sanitize metadata information

    Arguments:
    metadata_df - a pandas dataframe storing at least a column indicated by gene_col
    gene_col - the column to sanitize
    control_barcodes - a list of entries in the gene column to ignore sanitation for

    Returns:
    The same metadata dataframe as input with a sanitized column
    """
    genes = metadata_df.loc[:, gene_col].squeeze()
    if len(genes) > 0:
        genes = [x.split("_")[0] if x not in control_barcodes else x for x in genes]

        metadata_df.loc[:, gene_col] = genes
        return metadata_df.copy()
    else:
        return metadata_df
