import pandas as pd
from pycytominer import aggregate
from pycytominer.cyto_utils import infer_cp_features


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
    genes = [x.split("_")[0] if x not in control_barcodes else x for x in genes]

    metadata_df.loc[:, gene_col] = genes
    return metadata_df.copy()


def aggregate_pooled(site_df, strata, features, operation):
    """
    Aggregate each pooled cell painting site independently

    The arguments are the same as provided to pycytominer.aggregate()

    Returns:
    A dataframe of aggregated site-level profiles with perturbation count
    """

    site_counts_df = (
        site_df.loc[:, strata]
        .value_counts()
        .reset_index()
        .rename(columns={0: "Metadata_pert_count"})
    )

    site_df = aggregate(site_df, strata=strata, features=features, operation=operation)

    site_df = site_counts_df.merge(site_df, left_on=strata, right_on=strata)

    return site_df


def get_prop(df_group, pert_col="Metadata_pert_count"):
    """
    Get the proportion of total single cells that specific site contributed.
    Applied to a pandas.DataFrame.groupby() on the aggregated column level
    """
    pert_count = df_group[pert_col]
    prop_result = pd.DataFrame(pert_count / pert_count.sum()).rename(
        columns={pert_col: "Metadata_pert_prop"}
    )
    return prop_result.merge(df_group, left_index=True, right_index=True)


def approx_aggregate_piecewise(df, agg_cols, pert_col="Metadata_pert_count"):
    """
    Given previously aggregated files concatenated with a column indicating count,
    further aggregate.
    """
    meta_features = infer_cp_features(df, metadata=True)
    cp_features = df.drop(meta_features, axis="columns").columns.tolist()

    agg_df = (
        df.groupby(agg_cols)
        .apply(lambda x: get_prop(x, pert_col=pert_col))
        .reset_index(drop=True)
    )

    agg_df = agg_df.loc[:, meta_features].merge(
        agg_df.loc[:, cp_features].multiply(agg_df.Metadata_pert_prop, axis="rows"),
        left_index=True,
        right_index=True,
    )

    agg_df = agg_df.groupby(agg_cols).sum().reset_index()
    return agg_df
