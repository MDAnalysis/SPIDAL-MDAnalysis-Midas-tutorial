# combine output into full symmetric distance matrix

import json
import numpy as np

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

    print("Missing data for: {}".format(missing))

    return D

if __name__ == "__main__":

    outfilename = "distance_matrix.npz"
    traj = "trajectories.json"
    manifest = "manifest.json"

    D = combine(traj, manifest)

    np.save(outfilename, D)
    print("Created full distance matrix '{0}'".format(outfilename))


