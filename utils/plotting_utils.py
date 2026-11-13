import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn.objects as so
from matplotlib.patches import Circle, RegularPolygon
from matplotlib.path import Path
from matplotlib.projections import register_projection
from matplotlib.projections.polar import PolarAxes
from matplotlib.spines import Spine
from matplotlib.transforms import Affine2D


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

def radar_factory(num_vars, frame='circle'):
    """
    Create a radar chart with `num_vars` Axes.

    This function creates a RadarAxes projection and registers it.

    Parameters
    ----------
    num_vars : int
        Number of variables for radar chart.
    frame : {'circle', 'polygon'}
        Shape of frame surrounding Axes.

    """
    # calculate evenly-spaced axis angles
    theta = np.linspace(0, 2*np.pi, num_vars, endpoint=False)

    class RadarTransform(PolarAxes.PolarTransform):

        def transform_path_non_affine(self, path):
            # Paths with non-unit interpolation steps correspond to gridlines,
            # in which case we force interpolation (to defeat PolarTransform's
            # autoconversion to circular arcs).
            if path._interpolation_steps > 1:
                path = path.interpolated(num_vars)
            return Path(self.transform(path.vertices), path.codes)

    class RadarAxes(PolarAxes):

        name = f'radar{num_vars}'
        PolarTransform = RadarTransform

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            # rotate plot such that the first axis is at the top
            self.set_theta_zero_location('N')

        def fill(self, *args, closed=True, **kwargs):
            """Override fill so that line is closed by default"""
            return super().fill(closed=closed, *args, **kwargs)

        def plot(self, *args, **kwargs):
            """Override plot so that line is closed by default"""
            lines = super().plot(*args, **kwargs)
            for line in lines:
                self._close_line(line)

        def _close_line(self, line):
            x, y = line.get_data()
            # FIXME: markers at x[0], y[0] get doubled-up
            if x[0] != x[-1]:
                x = np.append(x, x[0])
                y = np.append(y, y[0])
                line.set_data(x, y)

        def set_varlabels(self, labels):
            self.set_thetagrids(np.degrees(theta), labels)

        def _gen_axes_patch(self):
            # The Axes patch must be centered at (0.5, 0.5) and of radius 0.5
            # in axes coordinates.
            if frame == 'circle':
                return Circle((0.5, 0.5), 0.5)
            elif frame == 'polygon':
                return RegularPolygon((0.5, 0.5), num_vars,
                                      radius=.5, edgecolor="k")
            else:
                raise ValueError("Unknown value for 'frame': %s" % frame)

        def _gen_axes_spines(self):
            if frame == 'circle':
                return super()._gen_axes_spines()
            elif frame == 'polygon':
                # spine_type must be 'left'/'right'/'top'/'bottom'/'circle'.
                spine = Spine(axes=self,
                              spine_type='circle',
                              path=Path.unit_regular_polygon(num_vars))
                # unit_regular_polygon gives a polygon of radius 1 centered at
                # (0, 0) but we want a polygon of radius 0.5 centered at (0.5,
                # 0.5) in axes coordinates.
                spine.set_transform(Affine2D().scale(.5).translate(.5, .5)
                                    + self.transAxes)
                return {'polar': spine}
            else:
                raise ValueError("Unknown value for 'frame': %s" % frame)

    register_projection(RadarAxes)
    return theta

def update_projection(axes, ax, projection='3d', fig=None, varlabels=None):
    if fig is None:
        fig = plt.gcf()
    rows, cols, start, stop = ax.get_subplotspec().get_geometry()
    axes.flat[start].remove()
    axes.flat[start] = fig.add_subplot(rows, cols, start+1, projection=projection)
    if varlabels:
        axes.flat[start].set_varlabels(varlabels)

def radar_plot(subsets,factory_names,factories,df,file_name):
    fig, axs = plt.subplots(figsize=(9, 9), nrows=1, ncols=3)
    for ax, subset, sub_factory in zip(axs.flat,subsets.keys(), factory_names):
        ax.set_title(subset)
        update_projection(axs,ax,sub_factory,fig,subsets[subset])
    for ax, subset,factory in zip(axs.flat,subsets.keys(),factories):
        fns = [np.mean,lambda x:np.percentile(x,10),lambda x:np.percentile(x,25),np.median,lambda x:np.percentile(x,75),lambda x:np.percentile(x,90)]
        for fn in fns:
            sub_data = []
            for each_subset in subsets[subset]:
                sub_data.append(fn(df[[x for x in df.columns if each_subset in x]].values))
            ax.plot(factory,sub_data)

    axs.flat[2].legend(['Mean', '10th percentile','Lower Quartile','Median','Upper Quartile', '90th percentile'],loc=(0.9, .95),
                            labelspacing=0.1, fontsize='small')
    fig.savefig(file_name, dpi=300, bbox_inches="tight", pad_inches=0.3)