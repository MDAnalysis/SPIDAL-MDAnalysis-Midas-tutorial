#!/usr/bin/env python
# Perform Path Similarity Analysis on a submatrix

"""Perform Path Similarity Analysis (PSA) using MDAnalysis.analysis.psa.

Provide all trajectories and topology files in a JSON input file (two
lists, one with topology files, the other with trajectory files) or on
the command line. There *must* be one topology file TOPOL for each
trajectory file TRAJ. The `--nsplit N` argument is required: it
indicates where to split the list of trajectories T so that the
distance matrix between the trajectories in `T[:N]` and `T[N:]`.

"""


import sys
import time
import argparse
from collections import OrderedDict
import itertools

import json
import numpy as np

import MDAnalysis as mda
from MDAnalysis.analysis import psa

def psa_partial(universesA, universesB, metric="discrete_frechet", selection="name CA"):
    """Calculate path similarity metric between lists of universes.

    Arguments
    ---------
    universesA, universesB : lists
         lists of MDAnalysis.Universe instances
    metric : string, optional
         label of the metric or distance function, can be one of "discrete_frechet",
         "hausdorff", "weighted_average_hausdorff", "average_hausdorff",
         "hausdorff_neighbors". Note that only "discrete_frechet" and "hausdorff"
         are true metrics.
    selection : string, optional
         MDAnalysis selection string that, when applied to *all* universes generates
         the subset of atoms that should be compared.

    Returns
    -------
    numpy.array(len(universesA), len(universesB)) : distance matrix
         The matrix of all the mutual distances D[i, j] with 0 <= i <
         len(A) and 0 <= j < len(B)

    Note
    ----
    Each universe is transferred to memory
    (`Universe.transfer_to_memory()`) and this permanently changes the
    universe, even in the calling code. The advantage is that
    subsequent calls to `transfer_to_memory()` are no-ops. However,
    for very big trajectories, memory problems might occur.

    """
    _metric = psa.get_path_metric_func(metric)

    # transfer to memory for instantaneous time series extraction
    # (typically not a problem for these trajectories and submatrix sizes; if memory becomes an
    # issue then do not transfer to memory and store intermediate trajectories on disk)
    for u in itertools.chain(universesA, universesB):
        u.transfer_to_memory()

    # submatrix of d[i, j] with i from A and j from B
    D = np.zeros((len(universesA), len(universesB)))

    # not optimized: could transpose to keep larger axis outside,
    # cache results (compare u_i and u_j), and generate 0 for u_i == u_j
    for i, u_i in enumerate(universesA):
        g1 = u_i.select_atoms(selection)
        P = u_i.trajectory.timeseries(asel=g1, format="fac")
        for j, u_j in enumerate(universesB):
            g2 = u_j.select_atoms(selection)
            Q = u_j.trajectory.timeseries(asel=g2, format="fac")

            # compute distance between paths
            D[i, j] = _metric(P, Q)

    return D

class StopWatch(OrderedDict):
    fmt = "{0:20s}  {1:8.3f} s"

    def tic(self, label):
        if label in self:
            raise ValueError("label {} already exists".format(label))
        self[label] = time.time()

    def show(self):
        print("----------[ TIMING ]--------------------")
        labels = self.keys()
        start = labels[0]
        for stop in labels[1:]:
            dt = self[stop] - self[start]
            print(self.fmt.format(stop, dt))
            start = stop
        print("----------------------------------------")
        print(self.fmt.format("total time",
                              self[labels[-1]] - self[labels[0]]))
        print("----------------------------------------")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-n", "--nsplit", type=int, required=True, default=None,
                        metavar="N",
                        help="split list of trajectories T so that the distance "
                        "submatrix D[i, j] will be calculated over the cartesian "
                        "product T[:N] x T[N:].")
    parser.add_argument("-f", "--inputfile",
                        help="JSON file with lists of topologies and "
                        "trajectories (or use --trajectories/--topologies)")
    parser.add_argument("--topologies", required=False, nargs="+",
                        metavar="TOPOLOGY",
                        help="List of topology files"),
    parser.add_argument("--trajectories", required=False, nargs="+",
                        metavar="TRAJ",
                        help="List of trajectory files")
    parser.add_argument("-o", "--outfile", default="distances.npy",
                        help="Distance matrix in numpy format")
    args = parser.parse_args()

    if args.inputfile:
        print("Loading paths from JSON file {}".format(args.inputfile))
        with open(args.inputfile) as inp:
            topologies, trajectories = json.load(inp)
    else:
        print("Using paths from command line")
        topologies = args.topologies
        trajectories = args.trajectories

    if len(topologies) != len(trajectories):
        raise ValueError("Need exactly one topology file for each trajectory")

    print("Processing {} trajectories.".format(len(trajectories)))
    print("Splitting trajectories in two blocks of length {0} and {1}".format(
         args.nsplit, len(trajectories) - args.nsplit))

    timer = StopWatch()
    timer.tic('init')

    # load trajectories
    universes = [mda.Universe(topology, trajectory) for topology, trajectory in
                 itertools.izip(topologies, trajectories)]
    timer.tic("load Universes")

    # run distance calculation and produce submatrix
         
    uA = universes[:args.nsplit]
    uB = universes[args.nsplit:]
    print("Calculating D (shape {0} x {1}) with {2} entries".format(
         len(uA), len(uB), len(uA) * len(uB)))

    D = psa_partial(uA, uB, metric="discrete_frechet", selection="name CA")

    timer.tic("PSA distance matrix")

    np.save(args.outfile, D)
    timer.tic("saving output")

    timer.show()


