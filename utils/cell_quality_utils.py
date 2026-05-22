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
1 - Perfect Cell (same as Simple)
2 - Great Cell (same as Simple)
3 - Imperfect High
  * This cell has many different barcodes mapped with various scores
  * The top scoring barcode has a 100% score and also has the most number of barcodes
4 - Imperfect Low
  * This cell has many different barcodes mapped with various scores
  * The top scoring barcode has <100% score and also has the most number of barcodes
5 - Bad Cell (same as Simple)

"""

import pandas as pd
import numpy as np
import sys
from functools import reduce


class CellQuality:
    def __init__(
        self,
        method,
        avg_col="mean",
        count_col="count",
        category_col_index="Metadata_Quality_Index",
        category_class_name="Metadata_Quality_Name",
    ):
        self.method = method
        self.avg_col = avg_col
        self.count_col = count_col
        self.category_col_index = category_col_index
        self.category_class_name = category_class_name

        if self.method == "simple":
            self.categorize = simple_categorize
        elif self.method == "simple_T7":
            self.categorize = simple_T7_categorize
        elif self.method == "simple_plus":
            self.categorize = simple_plus_categorize
        elif self.method == "multiples_T7":
            self.categorize = multiples_T7_categorize

        category_dict = get_cell_quality_dict(self.method)
        self.category_df = (
            pd.DataFrame(category_dict, index=[self.category_class_name])
            .transpose()
            .reset_index()
            .rename({"index": self.category_col_index}, axis="columns")
        )

    def assign_cell_quality(self, assigned_spot_df, parent_col, score_col, gene_col, barcode_col, SBScycles=12, match_to_library=False, library_structure=None):
        if len(library_structure) > 1:
            crispr_barcode_gene_dict = {}
            for col_to_match in library_structure.keys():
                barcode_group = assigned_spot_df.groupby(
                    [parent_col] + [f"{gene_col}_{col_to_match}", f"{barcode_col}_{col_to_match}"]
                )[f"{score_col}_{col_to_match}"]
                mergedf = pd.merge(
                    barcode_group.mean().reset_index(),
                    barcode_group.count().reset_index(),
                    on=[parent_col] + [f"{gene_col}_{col_to_match}", f"{barcode_col}_{col_to_match}"],
                    suffixes=["_mean", "_count"],
                )
                crispr_barcode_gene_dict[col_to_match] = mergedf
            crispr_barcode_gene_df = reduce(lambda left, right: pd.merge(left, right, how='outer'), crispr_barcode_gene_dict.values())
            
            cols_to_match = [f"{gene_col}_{col_to_match}" for col_to_match in library_structure.keys()]
            # Note if the genes matched to each iBAR match
            assigned_spot_df[f'GeneCallsMatch'] = assigned_spot_df[cols_to_match].eq(assigned_spot_df[cols_to_match[0]], axis=0).all(axis=1)
            # Assign a consensus barcode and gene if iBARs match
            assigned_spot_df['Barcode_MatchedTo_Barcode'] = "Unmatched" #TODO abstract col name
            assigned_spot_df['Barcode_MatchedTo_GeneCode'] = "Unmatched"
            assigned_spot_df.loc[assigned_spot_df['GeneCallsMatch'] == True, "Barcode_MatchedTo_Barcode"] = assigned_spot_df.loc[assigned_spot_df['GeneCallsMatch'] == True, cols_to_match[0]]
            assigned_spot_df.loc[assigned_spot_df['GeneCallsMatch'] == True, "Barcode_MatchedTo_GeneCode"] = assigned_spot_df.loc[assigned_spot_df['GeneCallsMatch'] == True, cols_to_match[0].replace("GeneCode","Barcode")]
            quality_df = (
                pd.DataFrame(
                    assigned_spot_df.groupby(parent_col).apply(
                        lambda x: self.categorize(x, cols_to_match[1]),
                        include_groups=False,
                    ),
                    columns=[self.category_col_index],
                )
                .reset_index()
                .merge(crispr_barcode_gene_df, on=parent_col)
                .merge(assigned_spot_df[[parent_col, "Barcode_MatchedTo_Barcode", "Barcode_MatchedTo_GeneCode"]], on=parent_col)
            ).assign(Quality_Method=self.method)
        else:
            barcode_group = assigned_spot_df.groupby(
                [parent_col] + [gene_col, barcode_col]
            )[score_col]

            crispr_barcode_gene_df = pd.merge(
                barcode_group.mean().reset_index(),
                barcode_group.count().reset_index(),
                on=[parent_col] + [gene_col, barcode_col],
                suffixes=["_mean", "_count"],
            )
            quality_df = (
                pd.DataFrame(
                    crispr_barcode_gene_df.groupby(parent_col).apply(
                        lambda x: self.categorize(
                            x, score_col=score_col, SBScycles=SBScycles
                        ),
                        include_groups=False,
                    ),
                    columns=[self.category_col_index],
                )
                .reset_index()
                .merge(crispr_barcode_gene_df, on=parent_col)
            ).assign(Quality_Method=self.method)

        cell_quality_dict = get_cell_quality_dict(self.method)
        quality_df[self.category_class_name] = quality_df[self.category_col_index].map(
            cell_quality_dict
        )
        return quality_df

    def summarize_cell_quality_counts(self, quality_df, parent_col):
        dup_cols = [parent_col, self.category_col_index]
        quality_count_df = (
            quality_df.drop_duplicates(subset=dup_cols)
            .loc[:, self.category_col_index]
            .value_counts()
            .reset_index()
            .rename(
                {
                    "count": "Cell_Count",
                },
                axis="columns",
            )
            .merge(self.category_df, on=self.category_col_index)
        )
        return quality_count_df

    def summarize_perturbation_quality_counts(
        self, quality_df, parent_col, group_cols, guide=False
    ):

        category_group_cols = group_cols + [self.category_col_index]
        category_group_cols = list(set(category_group_cols))

        if guide:
            level = "Guide"
        else:
            level = "Gene"

        summary_df = (
            quality_df.groupby(category_group_cols)[parent_col]
            .count()
            .reset_index()
            .rename({parent_col: f"Cell_Count_Per_{level}"}, axis="columns")
            .merge(self.category_df, on=self.category_col_index, how="left")
        )
        return summary_df


def get_cell_quality_dict(method):
    cell_quality_dict = {
        "simple": {1: "Perfect", 2: "Great", 3: "Imperfect", 4: "Bad"},
        "simple_T7": {
            1: "Perfect",
            2: "Great",
            3: "Imperfect_LowConfidence",
            4: "Imperfect_Multiples",
            5: "Bad",
        },
        "simple_plus": {
            1: "Perfect",
            2: "Great",
            3: "Imperfect-High",
            4: "Imperfect-Low",
            5: "Bad",
        },
        "multiples_T7": {
            1: "Perfect",
            2: "Good",
            3: "Acceptable",
            4: "Bad"
        },
    }
    return cell_quality_dict[method]


def filter_to_top_BC(SBS_df, parent_compartment, SBS_score_col_mean):
    SBS_df = SBS_df.loc[
        SBS_df.groupby(f"Parent_{parent_compartment}")[SBS_score_col_mean].idxmax()
    ].reset_index(drop=True)
    return SBS_df


def simple_categorize(
    parent_cell, score_col, avg_col="mean", count_col="count", SBScycles=12
):
    if len(score_col) > 1:
        print(f"simple_categorize method incompatible with matching multiple barcodes")
        sys.exit(1)
    # Written for SBS methods that produce many SBS foci/cell
    score_col_avg = f"{score_col}_{avg_col}"
    count_col_avg = f"{score_col}_{count_col}"

    parent_cell = parent_cell.sort_values(score_col_avg, ascending=False).reset_index(
        drop=True
    )

    num_barcodes = parent_cell.shape[0]
    max_score = max(parent_cell[score_col_avg])
    max_count = max(parent_cell[count_col_avg])

    # Only one barcode identified in the cell (any number of times)
    if num_barcodes == 1:
        # If all the calls in the cell were perfect calls then Cell is "Perfect"
        if max_score == 1:
            score = 1
        # If the calls in the cell were not all perfect calls then Cell is "Great"
        else:
            score = 2
    # Multiple different barcodes are read in the same cell
    else:
        list_bc_with_max_score = parent_cell.index[
            parent_cell[score_col_avg] == max_score
        ].values
        list_bc_most_identified = parent_cell.index[
            parent_cell[count_col_avg] == max_count
        ].values

        # If the top score is shared by multiple barcodes then Cell is "Bad"
        if len(list_bc_with_max_score) != 1:
            score = 4
        # If only one barcode is top scoring
        else:
            # If there is not a single barcode that is most identified
            if len(list_bc_most_identified) != 1:
                score = 4
            else:
                # If the top scoring barcode is the most common barcode
                if list_bc_with_max_score[0] == list_bc_most_identified[0]:
                    score = 3
                # If the top scoring barcode is not the most common barcode
                else:
                    score = 4
    return score


def simple_T7_categorize(parent_cell, score_col, avg_col="mean", SBScycles=12):
    if len(score_col) > 1:
        print(f"simple_T7_categorize method incompatible with matching multiple barcodes")
        sys.exit(1)
    # Written for SBS methods that produce median of 1 SBS focus/cell
    score_col_avg = f"{score_col}_{avg_col}"
    parent_cell = parent_cell.sort_values(score_col_avg, ascending=False).reset_index(
        drop=True
    )
    num_barcodes = parent_cell.shape[0]
    max_score = max(parent_cell[score_col_avg])
    ham1score = (SBScycles - 1) / (SBScycles)

    # Only one barcode identified in the cell (any number of times)
    if num_barcodes == 1:
        # If the call/s in the cell were perfect then Cell is "Perfect"
        if max_score == 1:
            score = 1
        # If the call/s in the cell have average score of off-by-one or better then Cell is "Great"
        elif ham1score >= max_score < 1:
            score = 2
        # If the call/s are on average worse than off-by-one then the Cell is "Imperfect_LowConfidence"
        else:
            score = 3
    # Multiple different barcodes are read in the same cell
    else:
        # If there is a Perfect call and the other/s is low confidence then Cell is "Imperfect_Multiples"
        if max_score == 1:
            other_scores = [x for x in parent_cell[score_col_avg] if x != 1]
            if all([x < ham1score for x in other_scores]):
                score = 4
            else:
                score = 5
        else:
            score = 5
    return score


def simple_plus_categorize(
    parent_cell, score_col, avg_col="mean", count_col="count", SBScycles=12
):
    if len(score_col) > 1:
        print(f"simple_plus_categorize method incompatible with matching multiple barcodes")
        sys.exit(1)
    score_col_avg = f"{score_col}_{avg_col}"
    count_col_avg = f"{score_col}_{count_col}"

    parent_cell = parent_cell.sort_values(score_col_avg, ascending=False).reset_index(
        drop=True
    )

    num_barcodes = parent_cell.shape[0]
    max_score = max(parent_cell[score_col_avg])
    max_count = max(parent_cell[count_col_avg])

    if num_barcodes == 1:
        if max_score == 1:
            score = 1
        else:
            score = 2
    else:
        list_bc_with_max_score = parent_cell.index[
            parent_cell[score_col_avg] == max_score
        ].values
        list_bc_most_identified = parent_cell.index[
            parent_cell[count_col_avg] == max_count
        ].values

        if len(list_bc_with_max_score) != 1:
            score = 5
        else:
            if len(list_bc_most_identified) != 1:
                score = 5
            else:
                if list_bc_with_max_score[0] == list_bc_most_identified[0]:
                    if max_score == 1:
                        score = 3
                    else:
                        score = 4
                else:
                    score = 5
    return score

def multiples_T7_categorize(parent_cell,gene_col1):
    # Written for SBS methods that produce median of 1 SBS focus/cell AND have multiple barcodes (iBARs)

    # TODO this should be abstracted more
    if parent_cell[gene_col1].nunique() == 1:
        if all(parent_cell['Spot_Category'] == "Perfect"):
            score = 1
        elif all(parent_cell['Spot_Category'].isin(["Perfect", "Good"])):
            score = 2
        elif all(parent_cell['Spot_Category'].isin(["Perfect", "Good", "Acceptable"])):
            score = 3
        # any bad spots, any recombinant spots
        else:
            score = 4
    else:
        score = 4
    return score