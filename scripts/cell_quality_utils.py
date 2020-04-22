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
    def __init__(self, method, avg_col="mean", count_col="count"):
        self.method = method
        self.avg_col = avg_col
        self.count_col = count_col

        if self.method == "simple":
            self.categorize = simple_categorize
        elif self.method == "simple_plus":
            self.categorize = simple_plus_categorize

    def define_cell_quality(self):
        return get_cell_quality_dict(method=self.method)


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
    }

    return cell_quality_dict[method]


def simple_categorize(parent_cell, score_col, avg_col="mean", count_col="count"):

    score_col_avg = f"{score_col}_{avg_col}"
    count_col_avg = f"{score_col}_{count_col}"

    parent_cell = parent_cell.sort_values(score_col_avg, ascending=False).reset_index(
        drop=True
    )

    num_barcodes = parent_cell.shape[0]
    max_score = max(parent_cell.Barcode_MatchedTo_Score_mean)
    max_count = max(parent_cell.Barcode_MatchedTo_Score_count)

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


def simple_plus_categorize(parent_cell, score_col, avg_col="mean", count_col="count"):

    score_col_avg = f"{score_col}_{avg_col}"
    count_col_avg = f"{score_col}_{count_col}"

    parent_cell = parent_cell.sort_values(score_col_avg, ascending=False).reset_index(
        drop=True
    )

    num_barcodes = parent_cell.shape[0]
    max_score = max(parent_cell.Barcode_MatchedTo_Score_mean)
    max_count = max(parent_cell.Barcode_MatchedTo_Score_count)

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
