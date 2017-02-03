#!/usr/bin/env python
"""
Generate lists of input files for the example data.
"""
import sys
import os
import glob
import json

DATADIR = "tutorial-data"

DIRECTORIES = {
    'DIMS':  os.path.join(DATADIR, "DIMS"),
    'FRODA': os.path.join(DATADIR, "FRODA"),
}



def make_file_list(method, topdir, start=None, stop=None, step=None):
    """Lists of corresponding topology and trajectory files.

    Match trajectories to appropriate topology file.
    """

    # contains hard coded paths

    trj_paths = glob.glob(
        os.path.join(topdir, DIRECTORIES[method], "trajectories", "*.dcd"))
    top_paths = {
        'DIMS': os.path.join(topdir, DIRECTORIES['DIMS'], "topologies", "adk4ake.psf"),
        'FRODA': os.path.join(topdir, DIRECTORIES['FRODA'], "topologies", "1ake.pdb")
    }

    trj = trj_paths[start:stop:step]
    top = len(trj) * [top_paths[method]]

    if len(trj) == 0:
	raise OSError("No trajectories for {0} in {1}".format(method, os.path.join(topdir, DIRECTORIES[method], "trajectories")))
    
    if not os.path.exists(top_paths[method]):
	 raise OSError("No topology for {0} in {1}".format(method, top_paths[method]))

    return top, trj

def echo(*args):
    sys.stderr.write("".join(args) + "\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),
                        default=sys.stdout,
                        help="write topology and trajectory file lists in JSON to "
                        "OUTFILE or stdout if no filename was supplied")
    parser.add_argument("-T", "--topdir", default=os.curdir,
                        help="Path to the directory containing '{}'".format(DATADIR))
    parser.add_argument("-b", "--start", type=int, default=None,
                        help="index of first trajectory, default  None, i.e., 0")
    parser.add_argument("-e", "--stop", type=int, default=None,
                                        help="index of last trajectory, default None, i.e., last")
    parser.add_argument("-s", "--step", type=int, default=None,
            help="step across start:stop:step, default None, i.e., 1")
    args = parser.parse_args()

    toppath = os.path.abspath(args.topdir)

    if not os.path.exists(toppath):
        raise ValueError("Path for --topdir {} not found".format(toppath))
    if not os.path.exists(os.path.join(toppath, DATADIR)):
        raise ValueError("Path for --topdir {0} does NOT contain {1}".format(toppath, DATADIR))
    echo("TOPPATH = {}".format(toppath))
    echo("DATADIR = {}".format(os.path.join(toppath, DATADIR)))

    files = make_file_list('DIMS', toppath,
                           start=args.start, stop=args.stop, step=args.step)
    f2 = make_file_list('FRODA', toppath,
                        start=args.start, stop=args.stop, step=args.step)
    files[0].extend(f2[0])
    files[1].extend(f2[1])
    del f2

    try:
        json.dump(files, args.outfile)
    finally:
        if args.outfile is not sys.stdout:
            args.outfile.close()
            echo("Wrote {0} files to {1}.".format(len(files[1]), args.outfile.name))


