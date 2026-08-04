import pandas as pd
import matplotlib.pyplot as plt
import seaborn.objects as so

def make_loc_df(sites_per_image_grid_side):
    # Creates x, y coordinates for plotting per-plate views.
    # Assumes image numbering starts in upper left corner and proceeds down
    final_order = []
    for i in range(1, sites_per_image_grid_side + 1):
        build_seq = list(
            zip(
                ([i] * (sites_per_image_grid_side + 1)),
                reversed(range(1, (sites_per_image_grid_side + 1))),
            )
        )
        final_order += build_seq
    sites_list = [
        *range(1, (sites_per_image_grid_side * sites_per_image_grid_side) + 1)
    ]
    loc_df = (
        pd.DataFrame(final_order)
        .rename(columns={0: "x_loc", 1: "y_loc"})
        .assign(Metadata_Site=sites_list)
        .astype({"Metadata_Site": int})
    )
    return loc_df

def apply_margin_titles(fig, col_order, row_order):
    # Replaces the "{col} | {row}" title seaborn.objects puts on every facet
    # with a single column label on the top row and a single row label on
    # the rightmost column, so each label appears once instead of once per subplot.
    numcols = len(col_order)
    numrows = len(row_order)
    for r in range(numrows):
        for c in range(numcols):
            ax = fig.axes[r * numcols + c]
            ax.set_title(str(col_order[c]) if r == 0 else "")
            if c == numcols - 1:
                ax.text(
                    1.05,
                    0.5,
                    str(row_order[r]),
                    transform=ax.transAxes,
                    rotation=270,
                    va="center",
                    ha="left",
                )

def make_plate_layout_plots(images_df, colorby, title, outpath, legend=False):
    numplates = images_df['Metadata_Plate'].nunique()
    numwells = images_df['Metadata_Well'].nunique()
    numsites = images_df['Metadata_Site'].nunique()
    pointsize = 40 if abs(numsites - 25) < abs(numsites - 100) else 20
    fig = plt.figure(figsize=(numplates*3,numwells*3))
    (
        so.Plot(images_df, x="x_loc", y="y_loc", color=colorby, text="Metadata_Site")
        .facet("Metadata_Plate","Metadata_Well")
        .add(so.Dot(pointsize=pointsize, marker="s"), legend=legend)
        .add(so.Text(color="w"))
        .label(x="", y="")
    ).on(fig).theme({"axes.facecolor": "w", "axes.edgecolor": "black"}).plot()
    for ax in fig.axes:
        ax.set_box_aspect(1)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_xticks([])
        ax.set_yticks([])
    if legend and fig.legends:
        old_legend = fig.legends[0]
        legend_title = old_legend.get_title().get_text()
        handles = old_legend.legend_handles
        labels = [text.get_text() for text in old_legend.get_texts()]
        fig.legends.remove(old_legend)
        for handle in handles:
            handle.set_sizes([400])
        # so.Plot's default legend is vertically centered on the whole figure,
        # which overlaps a subplot row whenever the grid has more than ~2 rows.
        # Anchoring to the top-right corner of the top-row axes keeps it level
        # with the top of the grid instead, clear of both the subplots and the title.
        fig.legend(
            handles,
            labels,
            loc="upper left",
            bbox_to_anchor=(1.05, 1.0),
            bbox_transform=fig.axes[0].transAxes,
        )
        title = f"{title}\n{legend_title}"
    fig.suptitle(title)
    # pad_inches guards against so.Plot's legend extent being computed slightly
    # smaller than its actual rendered size, which otherwise clips its right edge
    fig.savefig(outpath, dpi=300, bbox_inches="tight", pad_inches=0.3)
    plt.close(fig)