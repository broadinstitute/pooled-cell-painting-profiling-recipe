"""
Functions to determine cell quality

# Method: Simple

1 - Perfect Cell
  * This cell has a single or multiple of the same barcodes mapped, all with 100% score
2 - Great Cell
  * This cell has multiple of the same barcodes mapped
  * An average score less than 100%
3 - Imperfect Cell
  * This cell has many different barcodes mapped with various scores
  * The top scoring barcode also has the most number of barcodes
4 - Bad Cell
  * This cell has many different barcodes mapped with various scores
  * The top scoring barcode has fewer spots than any other barcode assigned to the cell

# Method: Simple Plus

This is the same as the "Simple" method, but splits Imperfect into two categories

3 - Imperfect High
  * This cell has many different barcodes mapped with various scores
  * The top scoring barcode has a 100% score and also has the most number of barcodes
4 - Imperfect Low
  * This cell has many different barcodes mapped with various scores
  * The top scoring barcode has <100% score and also has the most number of barcodes
5 - Bad Cell
  * Same as "Simple" method

"""

import pathlib
import pandas as pd


class CellQuality:
    def __init__(
        self,
        method,
        avg_col="mean",
        count_col="count",
        category_col_index="Metadata_Foci_Cell_Quality_Index",
        category_class_name="Metadata_Foci_Cell_Category",
    ):
        self.method = method
        self.avg_col = avg_col
        self.count_col = count_col
        self.category_col_index = category_col_index
        self.category_class_name = category_class_name

        if self.method == "simple":
            self.categorize = simple_categorize
        elif self.method == "simple_plus":
            self.categorize = simple_plus_categorize
        elif self.method == "feldman":
            self.categorize = feldman_categorize

        category_dict = self.define_cell_quality()
        self.category_df = (
            pd.DataFrame(category_dict, index=[self.category_class_name])
            .transpose()
            .reset_index()
            .rename({"index": self.category_col_index}, axis="columns")
        )

    def define_cell_quality(self):
        return get_cell_quality_dict(method=self.method)

    def assign_cell_quality(self, count_df, parent_cols, score_col, barcode_col):

        quality_estimate_df = (
            pd.DataFrame(
                count_df.groupby(parent_cols).apply(
                    lambda x: self.categorize(
                        x, score_col=score_col, barcode_col=barcode_col
                    )
                ),
                columns=[self.category_col_index],
            )
            .reset_index()
            .merge(count_df, on=parent_cols)
        ).assign(cell_quality_method=self.method)

        return quality_estimate_df

    def summarize_cell_quality_counts(self, quality_df, parent_cols):
        dup_cols = parent_cols + [self.category_col_index]

        quality_count_df = (
            quality_df.drop_duplicates(dup_cols)
            .loc[:, self.category_col_index]
            .value_counts()
            .reset_index()
            .rename(
                {
                    self.category_col_index: "Cell_Count",
                    "index": self.category_col_index,
                },
                axis="columns",
            )
            .merge(self.category_df, on=self.category_col_index)
        )

        return quality_count_df

    def summarize_perturbation_quality_counts(
        self, quality_df, parent_cols, group_cols, guide=False
    ):

        category_group_cols = group_cols + [self.category_col_index]
        category_group_cols = list(set(category_group_cols))

        if guide:
            level = "Guide"
        else:
            level = "Gene"

        summary_df = (
            quality_df.groupby(category_group_cols)[parent_cols[0]]
            .count()
            .reset_index()
            .rename({parent_cols[0]: f"Cell_Count_Per_{level}"}, axis="columns")
            .merge(self.category_df, on=self.category_col_index, how="left")
        )

        return summary_df


def get_cell_quality_dict(method):
    cell_quality_dict = {
        "simple": {1: "Perfect", 2: "Great", 3: "Imperfect", 4: "Bad"},
        "simple_plus": {
            1: "Perfect",
            2: "Great",
            3: "Imperfect-High",
            4: "Imperfect-Low",
            5: "Bad",
        },
        "feldman": {
            1: "Keep",
            2: "Toss_No_Perfect",
            3: "Toss_Multiple_Perfect",
            4: "Toss_Minority_Perfect",
        },
    }

    return cell_quality_dict[method]


def simple_categorize(
    parent_cell, score_col, barcode_col=None, avg_col="mean", count_col="count"
):

    score_col_avg = f"{score_col}_{avg_col}"
    count_col_avg = f"{score_col}_{count_col}"

    parent_cell = parent_cell.sort_values(score_col_avg, ascending=False).reset_index(
        drop=True
    )

    num_barcodes = parent_cell.shape[0]
    max_score = max(parent_cell.loc[:, score_col_avg])
    max_count = max(parent_cell.loc[:, count_col_avg])

    if num_barcodes == 1:
        if max_score == 1:
            score = 1
        else:
            score = 2
    else:
        max_score_idx = parent_cell.index[
            parent_cell[score_col_avg] == max_score
        ].values
        max_count_idx = parent_cell.index[
            parent_cell[count_col_avg] == max_count
        ].values

        if len(max_score_idx) != 1:
            score = 4
        else:
            if len(max_count_idx) != 1:
                score = 4
            else:
                if max_score_idx[0] == max_count_idx[0]:
                    score = 3
                else:
                    score = 4
    return score


def simple_plus_categorize(
    parent_cell, score_col, barcode_col=None, avg_col="mean", count_col="count"
):

    score_col_avg = f"{score_col}_{avg_col}"
    count_col_avg = f"{score_col}_{count_col}"

    parent_cell = parent_cell.sort_values(score_col_avg, ascending=False).reset_index(
        drop=True
    )

    num_barcodes = parent_cell.shape[0]
    max_score = max(parent_cell.loc[:, score_col_avg])
    max_count = max(parent_cell.loc[:, count_col_avg])

    if num_barcodes == 1:
        if max_score == 1:
            score = 1
        else:
            score = 2
    else:
        max_score_idx = parent_cell.index[
            parent_cell[score_col_avg] == max_score
        ].values
        max_count_idx = parent_cell.index[
            parent_cell[count_col_avg] == max_count
        ].values

        if len(max_score_idx) != 1:
            score = 5
        else:
            if len(max_count_idx) != 1:
                score = 5
            else:
                if max_score_idx[0] == max_count_idx[0]:
                    if max_score == 1:
                        score = 3
                    else:
                        score = 4
                else:
                    score = 5
    return score


def feldman_categorize(
    parent_cell,
    score_col,
    barcode_col="Barcode_MatchedTo_Barcode",
    avg_col=None,
    count_col=None,
):
    """
    Note that the Feldman categorize function works on non-averaged scores
    """
    num_barcodes = parent_cell.shape[0]
    max_score = max(parent_cell.loc[:, score_col])

    if max_score < 1:
        score = 2
    else:
        barcode_count_with_max_score = parent_cell.index[
            parent_cell[score_col] == max_score
        ].values

        perfect_matches = parent_cell.loc[barcode_count_with_max_score, :]

        if len(perfect_matches.loc[:, barcode_col].unique()) != 1:
            score = 3
        else:
            top_barcode_ratio = len(barcode_count_with_max_score) / num_barcodes
            if top_barcode_ratio > 0.5:
                score = 1
            else:
                score = 4

    return score
