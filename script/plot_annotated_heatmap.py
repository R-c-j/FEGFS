
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    font1 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': 24}
    font2 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': 24}
    ax.set_xticklabels(col_labels,font1)
    ax.set_yticklabels(row_labels,font2,rotation=30)

    # modify this part and change the label direction
    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=False, bottom=True,
                   labeltop=False, labelbottom=True)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=0.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar

def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=["black", "white"],
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A list or array of two color specifications.  The first is used for
        values below a threshold, the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts

# ARI
datasets = ["Pollen", "Goolam", "Patel", "Biase","Kelin","Zeisel"]
methods =["K-means\n+T-SNE","SC3","SSC","SIMLR","CORR","SinNLRR","Seurat","SSSC","FFMC"]
results = np.array([[0.839,0.932,0.838,0.938,0.845,0.932,0.850,0.826,0.955],
                    [0.520,0.687,0.544,0.490,0.508,0.910,0.883,0.544,0.910],
                    [0.511,0.789,0.706,0.757,0.774,0.656,0.888,0.412,0.877],
                    [0.720,0.922,1.000,0.740,1.000,0.883,0.581,1.000,0.860],
                    [0.581,0.798,0.049,0.639,0.743,0.627,0.819,0.033,0.781],
                    [0.525,0.807,0.035,0.704,0.672,0.805,0.516,0.359,0.781],
                    ])

# NMI
# datasets = ["Pollen", "Goolam", "Patel", "Biase","Kelin","Zeisel"]
# methods =["K-means\n+T-SNE","SC3","SSC","SIMLR","CORR","SinNLRR","Seurat","SSSC","FFMC"]
# results = np.array([[0.890,0.937,0.927,0.948,0.910,0.939,0.930,0.898,0.952],
#                     [0.681,0.844,0.720,0.650,0.664,0.881,0.872,0.706,0.881],
#                     [0.563,0.866,0.802,0.799,0.866,0.716,0.863,0.432,0.859],
#                     [0.767,0.917,1.000,0.821,1.000,0.848,0.690,1.000,0.878],
#                     [0.592,0.855,0.112,0.749,0.798,0.757,0.843,0.056,0.820],
#                     [0.669,0.754,0.119,0.765,0.699,0.766,0.700,0.422,0.774],
#                     ])

fig, ax = plt.subplots(figsize=(18,10))
im, cbar = heatmap(results, datasets, methods, ax=ax,
                   cmap=plt.get_cmap("GnBu"))
texts = annotate_heatmap(im, valfmt="{x:.3f}",size = 16)
font1 = {'family': 'Times New Roman',
         'weight': 'normal',
         'size': 24}
ax.set_title('ARI',fontsize=30,color='black',)
fig.tight_layout()

plt.savefig('./eps图/ARI.tif', dpi=300)
plt.show()
