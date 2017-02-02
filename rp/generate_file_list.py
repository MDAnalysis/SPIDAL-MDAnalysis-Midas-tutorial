#!/usr/bin/env python
"""
Generate lists of input files for the example data.
"""
import sys
import os
import glob
import json

DIRECTORIES = {
            "DIMS":
                "/nfs/homes3/sseyler/Projects/Enzymes/AdK/simulations/enhanced/DIMS/trj/1ake/core",
            "FRODA":
                "/nfs/homes3/sseyler/Simulations/adk/gp/trj/co/fit/core",
             }
trj_paths = dict((method, glob.glob(os.path.join(pth, "*.dcd"))) 
                         for method, pth in DIRECTORIES.items())
top_paths = {'DIMS': 
                 "/nfs/homes3/sseyler/Projects/Enzymes/AdK/simulations/enhanced/DIMS/top/adk4ake.psf",
             'FRODA': 
                 "/nfs/homes3/sseyler/Simulations/adk/gp/top/1ake.pdb"}

def make_file_list(method, start=None, stop=None, step=None):
    """Lists of corresponding topology and trajectory files"""

    trj = trj_paths[method][start:stop:step]
    top = len(trj) * [top_paths[method]]
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
    parser.add_argument("-b", "--start", type=int, default=None,
                        help="index of first trajectory, default  None, i.e., 0")
    parser.add_argument("-e", "--stop", type=int, default=None,
                                        help="index of last trajectory, default None, i.e., last")
    parser.add_argument("-s", "--step", type=int, default=None,
            help="step across start:stop:step, default None, i.e., 1")
    args = parser.parse_args()

    files = make_file_list('DIMS', start=args.start, stop=args.stop, step=args.step)
    f2 = make_file_list('FRODA', start=args.start, stop=args.stop, step=args.step)
    files[0].extend(f2[0])
    files[1].extend(f2[1])
    del f2

    try:
        json.dump(files, args.outfile)
    finally:
        if args.outfile is not sys.stdout:
            args.outfile.close()
            echo("Wrote {0} files to {1}.".format(len(files[1]), args.outfile.name))


