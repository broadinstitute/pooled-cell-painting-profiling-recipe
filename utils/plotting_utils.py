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

def make_plate_layout_plots(images_df, colorby, title, outpath, legend=False):
    numplates = images_df['Metadata_Plate'].nunique()
    numwells = images_df['Metadata_Well'].nunique()
    fig = plt.figure(figsize=(numwells*3,numplates*3))
    (
        so.Plot(images_df, x="x_loc", y="y_loc", color=colorby, text="Metadata_Site")
        .facet("Metadata_Well","Metadata_Plate")
        .add(so.Dot(pointsize=20, marker="s"), legend=legend)
        .add(so.Text(color="w"))
        .label(x="", y="")
    ).on(fig).theme({"axes.facecolor": "w", "axes.edgecolor": "black"}).plot()
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300,bbox_inches='tight')