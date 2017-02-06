#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""\
Combine output files from radical.pilot runs of PSA (block matrices)
into full symmetric distance matrix, save the complete matrix, and
cluster it (using Ward's linkage). The clustered matrix together with
its dendrogram is written to an output file.

If output files are still missing, warnings are printed and the
entries in the matrix are set to 0.

NOTE: This script is adapted to the tutorial and DIMS trajectories are
labelled in the plot with a thin bar whereas FRODA trajectories are
labelled with a heavy double bar. Change the code (look for the
function make_labels()) if you want to use it for other data.

"""

import os.path
import json
import numpy as np

import matplotlib.pyplot as plt

import matplotlib
import scipy.cluster.hierarchy


def combine(traj, manifest):
    trj = json.load(open(traj))
    N = len(trj[0])
    D = np.zeros((N, N))

    mf = json.load(open(manifest))

    missing = []
    for subD, _, (i0, i1), (j0, j1) in mf:
        try:
            D[i0:i1, j0:j1] = np.load(subD)
            D[j0:j1, i0:i1] = D[i0:i1, j0:j1].T
        except IOError:
            missing.append(subD)

    if missing:
        print("WARNING: Missing data: output files\n{}\n"
              "could not be found".format(", ".joind(missing)))

    return D

def cluster(distArray, method='ward', count_sort=False,
            distance_sort=False, no_plot=False, no_labels=True,
            color_threshold=4):
    """Cluster trajectories and optionally plot the dendrogram.

    :Arguments:
      *method*
        string, name of linkage criterion for clustering [``'ward'``]
      *no_plot*
        boolean, if ``True``, do not render the dendrogram [``False``]
      *no_labels*
        boolean, if ``True`` then do not label dendrogram [``True``]
      *color_threshold*
        For brevity, let t be the color_threshold. Colors all the
        descendent links below a cluster node k the same color if k is
        the first node below the cut threshold t. All links connecting
        nodes with distances greater than or equal to the threshold are
        colored blue. If t is less than or equal to zero, all nodes are
        colored blue. If color_threshold is None or ‘default’,
        corresponding with MATLAB(TM) behavior, the threshold is set to
        0.7*max(Z[:,2]). [``4``]]
    :Returns:
      list of indices representing the row-wise order of the objects
      after clustering
    """
    matplotlib.rcParams['lines.linewidth'] = 0.5

    Z = scipy.cluster.hierarchy.linkage(distArray, method=method)
    dgram = scipy.cluster.hierarchy.dendrogram(
        Z, no_labels=no_labels, orientation='left',
        count_sort=count_sort, distance_sort=distance_sort,
        no_plot=no_plot, color_threshold=color_threshold)
    return Z, dgram


def plot_clustered_distances(dist_matrix, filename=None, linkage='ward', count_sort=False,
                             distance_sort=False, labels=None,
                             figsize=4.5, labelsize=12):
    """Plot a clustered distance matrix.

    Using method *linkage* along with the corresponding
    dendrogram. Rows (and columns) are identified using the list of
    strings specified by `labels`.

    :Arguments:
      *filename*
         string, save figure to *filename* [``None``]
      *linkage*
         string, name of linkage criterion for clustering [``'ward'``]
      *count_sort*
        boolean, see scipy.cluster.hierarchy.dendrogram [``False``]
      *distance_sort*
        boolean, see scipy.cluster.hierarchy.dendrogram [``False``]
      *figsize*
         set the vertical size of plot in inches [``4.5``]
      *labels*
         list of labels, sorted like the original trajectories [``None``]
      *labelsize*
         set the font size for colorbar labels; font size for path labels on
         dendrogram default to 3 points smaller [``12``]

    If *filename* is supplied then the figure is also written to file (the
    suffix determines the file type, e.g. pdf, png, eps, ...). All other
    keyword arguments are passed on to :func:`pylab.imshow`.

    """

    # make sure that this is a square distance matrix
    assert dist_matrix.ndim == 2
    assert dist_matrix.shape[0] == dist_matrix.shape[1]

    npaths = len(dist_matrix)

    dgram_loc, hmap_loc, cbar_loc = _get_plot_obj_locs()
    aspect_ratio = 1.25

    fig = plt.figure(figsize=(figsize*aspect_ratio, figsize))

    ax_hmap = fig.add_axes(hmap_loc)
    ax_dgram = fig.add_axes(dgram_loc)

    Z, dgram = cluster(dist_matrix,
                       method=linkage,
                       count_sort=count_sort,
                       distance_sort=distance_sort)
    rowidx = colidx = dgram['leaves'] # get row-wise ordering from clustering
    ax_dgram.invert_yaxis() # Place origin at up left (from low left)

    minDist, maxDist = 0, np.max(dist_matrix)
    dist_matrix_clus = dist_matrix[rowidx,:]
    dist_matrix_clus = dist_matrix_clus[:,colidx]
    im = ax_hmap.matshow(dist_matrix_clus, aspect='auto', origin='lower',
                cmap=plt.cm.YlGn, vmin=minDist, vmax=maxDist)
    ax_hmap.invert_yaxis() # Place origin at upper left (from lower left)
    ax_hmap.locator_params(nbins=npaths)
    ax_hmap.set_xticks(np.arange(npaths), minor=True)
    ax_hmap.set_yticks(np.arange(npaths), minor=True)
    ax_hmap.tick_params(axis='x', which='both', labelleft='off',
                        labelright='off', labeltop='on', labelsize=0)
    ax_hmap.tick_params(axis='y', which='both', labelleft='on',
                        labelright='off', labeltop='off', labelsize=0)
    if labels:
        rowlabels = [labels[i] for i in rowidx]
        collabels = [labels[i] for i in colidx]
        ax_hmap.set_xticklabels(collabels, rotation='vertical',
                                size=(labelsize-4), multialignment='center', minor=True)
        ax_hmap.set_yticklabels(rowlabels, rotation='horizontal',
                                size=(labelsize-4), multialignment='left', ha='right',
                                minor=True)

    ax_color = fig.add_axes(cbar_loc)
    plt.colorbar(im, cax=ax_color, ticks=np.linspace(minDist, maxDist, 10),
                 format="%0.1f")
    ax_color.tick_params(labelsize=labelsize)

    # Remove major ticks from both heat map axes
    for tic in ax_hmap.xaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False
        tic.label1On = tic.label2On = False
    for tic in ax_hmap.yaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False
        tic.label1On = tic.label2On = False
    # Remove minor ticks from both heat map axes
    for tic in ax_hmap.xaxis.get_minor_ticks():
        tic.tick1On = tic.tick2On = False
    for tic in ax_hmap.yaxis.get_minor_ticks():
        tic.tick1On = tic.tick2On = False
    # Remove tickmarks from colorbar
    for tic in ax_color.yaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False

    if filename is not None:
        fig.savefig(filename, dpi=300, bbox_inches='tight')

    return Z, dgram, dist_matrix_clus

def _get_plot_obj_locs():
    """Find and return coordinates for dendrogram, heat map, and colorbar.

    :Returns:
      tuple of coordinates for placing the dendrogram, heat map, and
      colorbar in the plot.
    """
    plot_xstart = 0.04
    plot_ystart = 0.04
    label_margin = 0.155

    dgram_height = 0.2 # dendrogram heights(s)
    hmap_xstart = plot_xstart + dgram_height + label_margin

    # Set locations for dendrogram(s), matrix, and colorbar
    hmap_height = 0.8
    hmap_width = 0.6
    dgram_loc = [plot_xstart, plot_ystart, dgram_height, hmap_height]
    cbar_width = 0.02
    cbar_xstart = hmap_xstart + hmap_width + 0.01
    cbar_loc = [cbar_xstart, plot_ystart, cbar_width, hmap_height]
    hmap_loc =  [hmap_xstart, plot_ystart, hmap_width, hmap_height]

    return dgram_loc, hmap_loc, cbar_loc


def make_labels(traj, mapping, no_label=""):
    trajectories = json.load(open(traj))
    labels = []
    for trj in trajectories[1]:
        for pattern, label in mapping.items():
            if pattern in trj:
                labels.append(label)
                break
        else:
            labels.append(no_label)
    return labels


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--outfile", metavar="FILE",
                        default= "distance_matrix.npy",
                        help="save complete distance matrix to FILE")
    parser.add_argument("-p", "--plot", metavar="FILE",
                        default="psa_distance_matrix.pdf",
                        help="plot clustered distance matrix to FILE")
    parser.add_argument("-t", "--trajectories", metavar="FILE",
                        default="trajectories.json",
                        help="JSON file with the trajectories")
    parser.add_argument("-f", "--manifest", metavar="FILE",
                        default="manifest.json",
                        help="JSON file with the complete manifest of all files "
                        "including the indices of the block matrices")

    args = parser.parse_args()

    # reduce: combine all submatrices
    D = combine(args.trajectories, args.manifest)

    np.save(args.outfile, D)
    print("Created full distance matrix '{0}'".format(args.outfile))

    # analyze: plot clustered distance matrix
    # Label FRODA trajectories with a double bar and DIMS with a thin single bar
    plot_clustered_distances(D, filename=args.plot,
                             labels=make_labels(args.trajectories,
                                                {'DIMS': '.', 'FRODA': '--'}))
    print("Plotted clustered distance matrix '{0}'".format(args.plot))

    # plot as PNG for tutorial/web
    figpng = os.path.splitext(args.plot)[0] + ".png"
    plt.savefig(figpng, dpi=300)
    print("Plotted clustered distance matrix '{0}'".format(figpng))
