# radical.pilot script to run MDAnalysis Path Similarity Analysis
#
# =======================================
# configured for XSEDE TACC stampede
# =======================================
#
# USAGE: python rp_psa.py trajectories.json <blocksize> <cores> <session-name>
#
# ENVIRONMENT VARIABLES must be set:
#
#   export RADICAL_PILOT_PROJECT=TG-xxxxxx
#   export RADICAL_PILOT_DBURL="mongodb://user:pass@host:port/dbname"
#
#
# Provide a JSON file with two lists of files:
#
# [
#  [topol1, topol1, topol1, ..., topol2, ...],
#  [trj1, trj2, trj3, ...]
# ]
#
# where each trj file name has a corresponding topology filename;
# the topology filenames can (and should) be repeated as often as
# necessary.
#
# The entries must be absolute paths.
#
# For this example, no additional superposition is performed and
# the trajectories are used as is. CA atoms are selected.
#
# The Frechet distance matrix between ALL trajectories (all-vs-all)
# is calculated.

import os
from os.path import basename
import time

os.environ['RADICAL_PILOT_VERBOSE']='DEBUG'
# os.environ['RADICAL_PILOT_LOG_TGT']='/tmp/rp.log'

import sys
import radical.pilot as rp

import json

#----------------------------------------------------------------------------
if __name__ == "__main__":

    MY_STAGING_AREA = 'staging:///'
    SHARED_MDA_SCRIPT = 'mdanalysis_psa_partial.py'
    FILELIST = sys.argv[1]
    BLOCK_SIZE = int(sys.argv[2])
    cores = int(sys.argv[3])
    session_name = sys.argv[4]
    MANIFEST_NAME = "manifest.json"

    try:
        PROJECT = os.environ['RADICAL_PILOT_PROJECT']
        if not PROJECT:
            raise ValueError
    except (KeyError, ValueError):
        raise RuntimeError("Set RADICAL_PILOT_PROJECT to your XSEDE allocation")

    try:
        session   = rp.Session (name=session_name)
        c         = rp.Context ('ssh')
        session.add_context (c)

        print "initialize pilot manager ..."
        pmgr = rp.PilotManager (session=session)
        #pmgr.register_callback (pilot_state_cb)

        pdesc = rp.ComputePilotDescription ()
        pdesc.resource = "xsede.stampede"   # xsede.stampede or xsede.comet or ... check docs!
        pdesc.runtime  = 20 # minutes
        pdesc.cores    = cores
        pdesc.project  = PROJECT #Project allocation, pass through env var PROJECT
        pdesc.cleanup  = False
        pdesc.access_schema = 'ssh'

        # submit the pilot.
        pilot = pmgr.submit_pilots (pdesc)

        #initialize unit manager
        umgr = rp.UnitManager  (session=session, scheduler=rp.SCHED_DIRECT_SUBMISSION)

        #add pilot to unit manager
        umgr.add_pilots(pilot)

        # get ALL topology and trajectory files
        with open(FILELIST) as inp:
            topologies, trajectories = json.load(inp)

        # list of outputfiles <--> submatrix indices
        # (built during CU creation)
        manifest = []

        fshared_list   = []
        fname_stage = []
        # stage all files to the staging area
        src_url = 'file://%s/%s' % (os.getcwd(), SHARED_MDA_SCRIPT)

        #print src_url

        sd_pilot = {'source': src_url,
                    'target': os.path.join(MY_STAGING_AREA, SHARED_MDA_SCRIPT),
                    'action': rp.TRANSFER,
        }
        fname_stage.append(sd_pilot)

        # Synchronously stage the data to the pilot
        pilot.stage_in(fname_stage)

        # we create one CU for a block of the distance matrix
        cudesc_list = []

        for i in range(0, len(trajectories), BLOCK_SIZE):
            for j in range(i, len(trajectories), BLOCK_SIZE):
                fshared = list()
                shared = {'source': os.path.join(MY_STAGING_AREA, SHARED_MDA_SCRIPT),
                          'target': SHARED_MDA_SCRIPT,
                          'action': rp.LINK}
                fshared.append(shared)

                shared = [{'source': 'file://{0}'.format(trajectory),
                           'target' : basename(trajectory),
                           'action' : rp.LINK}
                              for trajectory in trajectories[i:i+BLOCK_SIZE]]
                fshared.extend(shared)
                if i!=j:
                    shared = [{'source': 'file://{0}'.format(trajectory),
                               'target' : basename(trajectory),
                               'action' : rp.LINK}
                                  for trajectory in trajectories[j:j+BLOCK_SIZE]]
                    fshared.extend(shared)
                # always copy all unique topology files
                shared = [{'source': 'file://{0}'.format(topology),
                           'target' : basename(topology),
                           'action' : rp.LINK}
                          for topology in set(topologies)]
                fshared.extend(shared)

                # block of topology / trajectory pairs
                #   block[:nsplit] + block[nsplit:]
                # The MDA script wants one long list of trajectories and the index nsplit
                # that indicates where to split the list to create the two groups of
                # trajectories that are compared against each other.
                block_top = topologies[i:i+BLOCK_SIZE] + topologies[j:j+BLOCK_SIZE]
                block_trj = trajectories[i:i+BLOCK_SIZE] + trajectories[j:j+BLOCK_SIZE]
                block = [block_top, block_trj]
                nsplit = len(trajectories[i:i+BLOCK_SIZE])
                imax = i + len(trajectories[i:i+BLOCK_SIZE])
                jmax = j + len(trajectories[j:j+BLOCK_SIZE])
                # should remember i, imax and j_jmax because we calculate the
                # submatrix D[i:i+di, j:j+dj] in this CU.
                block_json = "block-{0}-{1}__{2}-{3}.json".format(
                    i, imax, j, jmax)
                block_matrixfile = 'subdistances_{0}-{1}__{2}-{3}.npy'.format(
                    i, imax, j, jmax)

                manifest.append((block_matrixfile, block_json, (i, imax), (j, jmax)))

                # create input file for the cu and add share it
                with open(block_json, "w") as out:
                    json.dump(block, out)

                # share input json file
                shared = [{'source': 'file://{0}'.format(os.path.realpath(block_json)),
                           'target' : basename(block_json),
                           'action' : rp.LINK}
                ]
                fshared.extend(shared)

                # define the compute unit, to compute over the trajectory submatrix
                # TODO: store the offsets WITH the returned matrix (pkl or arrach archive) instead
                #       of encoding in filename
                cudesc = rp.ComputeUnitDescription()
                cudesc.executable     = "python"
                cudesc.pre_exec       = ["module load python; source activate mdaenv"] #Only for Stampede and with our conda env
                cudesc.input_staging  = fshared
                cudesc.output_staging = {'source': block_matrixfile,
                                         'target': block_matrixfile,
                                         'action': rp.TRANSFER}
                cudesc.arguments      = [SHARED_MDA_SCRIPT, '--nsplit', nsplit,
                                         '--inputfile', block_json,
                                         '--outfile', block_matrixfile, ]
                cudesc.cores          = 1

                cudesc_list.append(cudesc)

        # write manifest json file: use it later to piece submatrices back
        # together
        with open(MANIFEST_NAME, "w") as outfile:
            json.dump(manifest, outfile)
        print("Created manifest '{0}': (block_D, block_trj, (i, i+w), (j, j+w))".format(MANIFEST_NAME))

        # submit, run and wait and...
        #print "submit units to unit manager ..."
        units = umgr.submit_units(cudesc_list)

        #print "wait for units ..."
        umgr.wait_units()

        print "Creating Profile"
        with open('{0}.csv'.format(session.name),'w') as ProfFile:
            ProfFile.write('CU,New,StageIn,Allocate,Exec,StageOut,Done\n')
            for cu in units:
                timing_str=[cu.uid,'N/A','N/A','N/A','N/A','N/A','N/A']
                for states in cu.state_history:
                    if states.as_dict()['state']=='Scheduling':
                        timing_str[1]= (states.as_dict()['timestamp']-pilot.start_time).__str__()
                    elif states.as_dict()['state']=='AgentStagingInput':
                        timing_str[2]= (states.as_dict()['timestamp']-pilot.start_time).__str__()
                    elif states.as_dict()['state']=='Allocating':
                        timing_str[3]= (states.as_dict()['timestamp']-pilot.start_time).__str__()
                    elif states.as_dict()['state']=='Executing':
                        timing_str[4]= (states.as_dict()['timestamp']-pilot.start_time).__str__()
                    elif states.as_dict()['state']=='AgentStagingOutput':
                        timing_str[5]= (states.as_dict()['timestamp']-pilot.start_time).__str__()
                    elif states.as_dict()['state']=='Done':
                        timing_str[6]= (states.as_dict()['timestamp']-pilot.start_time).__str__()

                ProfFile.write(timing_str[0]+','+timing_str[1]+','+
                               timing_str[2]+','+timing_str[3]+','+
                               timing_str[4]+','+timing_str[5]+','+
                               timing_str[6]+'\n')


    except (KeyboardInterrupt, SystemExit) as e:
        # the callback called sys.exit(), and we can here catch the
        # corresponding KeyboardInterrupt exception for shutdown.  We also catch
        # SystemExit (which gets raised if the main threads exits for some other
        # reason).
        print "need to exit now: %s" % e

    except Exception as e:
        # Something unexpected happened in the pilot code above
        print "caught Exception: %s" % e
        raise

    finally :
        print "Closing session, exiting now ..."
        session.close(cleanup=False)

