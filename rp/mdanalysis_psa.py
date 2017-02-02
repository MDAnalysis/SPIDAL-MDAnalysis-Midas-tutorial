#!/usr/bin/env python
# Perform Path Similarity Analysis on a submatrix

"""Perform Path Similarity Analysis (PSA) using
MDAnalysis.analysis.psa.PSAnalysis. 

Provide all trajectories and topology files on the command line. There
*must* be one topology file TOPOL for each trajectory file TRAJ.

"""


import sys
import time
import argparse
from collections import OrderedDict

import json
import numpy as np

import MDAnalysis as mda
import MDAnalysis.analysis.psa

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
    parser.add_argument("--inputfile", 
                        help="JSON file with lists of topologies and "
                        "trajectories (or use --trajectories/--topologies)")
    parser.add_argument("--topologies", required=False, nargs="+",
                        metavar="TOPOLOGY",
                        help="List of topology files"), 
    parser.add_argument("--trajectories", required=False, nargs="+",
                        metavar="TRAJ",
                        help="List of trajectory files")
    parser.add_argument("--outfile", default="distances.npy",
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

    timer = StopWatch()
    timer.tic('init')
    
    # load trajectories
    universes = [mda.Universe(topology, trajectory) for topology, trajectory in
                 zip(topologies, trajectories)]    
    timer.tic("load Universes")

    # set up PSA
    P = MDAnalysis.analysis.psa.PSAnalysis(universes,
                                           ref_select="name CA", targetdir="tmp_psa")
    timer.tic("set up PSA")

    # generate paths (no superposition)
    P.generate_paths(align=False)
    timer.tic("generate PSA paths")

    # run distance calculation and produce submatrix
    P.run(metric="discrete_frechet")
    timer.tic("PSA distance matrix")
    
    np.save(args.outfile, P.D)
    timer.tic("saving output")    
              
    timer.show()

    
