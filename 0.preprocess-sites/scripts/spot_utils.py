import pathlib
import pandas as pd
import plotnine as gg
import matplotlib.pyplot as plt
import seaborn as sns


def spot_counts_per_cell_histogram(df, col, file, bins=50):
    plt.figure(num=None, figsize=(4, 3), dpi=300, facecolor="w", edgecolor="k")
    df.loc[:, col].squeeze().value_counts().hist(bins=bins)
    plt.xlabel("Number of Barcodes")
    plt.ylabel("Number of Cells")
    plt.tight_layout()
    plt.savefig(file)
    plt.close()


def spot_score_histogram(df, col, file, bins=50):
    plt.figure(num=None, figsize=(4, 3), dpi=300, facecolor="w", edgecolor="k")
    df.loc[:, col].squeeze().hist(bins=bins)
    plt.xlabel("Spot Scores (Alignment)")
    plt.ylabel("Number of Spots")
    plt.tight_layout()
    plt.savefig(file)
    plt.close()


def spot_count_score_jointplot(df, parent_col, score_col, file):
    avg_df = (
        pd.DataFrame(df.groupby(parent_col)[score_col].mean())
        .reset_index()
        .merge(
            (
                pd.DataFrame(df.loc[:, parent_col].squeeze().value_counts())
                .reset_index()
                .rename(
                    {parent_col: "Barcode_Count", "index": parent_col}, axis="columns"
                )
            ),
            on=parent_col,
        )
    )

    sns.jointplot(
        x=score_col,
        y="Barcode_Count",
        kind="reg",
        data=avg_df,
        scatter_kws={"s": 0.25},
        height=3,
    )
    plt.xlabel("Mean Spot Score per Cell (Alignment)")
    plt.ylabel("Number of Barcodes per Cell")
    plt.savefig(file, dpi=300)
    plt.close()


def category_counts(df, gene_cols, barcode_cols, score_cols, parent_cols, guide=False):
    ids = parent_cols + gene_cols
    if guide:
        ids += barcode_cols

    barcode_group = df.groupby(ids)[score_cols]

    count_df = pd.merge(
        barcode_group.mean().reset_index(),
        barcode_group.count().reset_index(),
        on=ids,
        suffixes=["_mean", "_count"],
    )

    return count_df
