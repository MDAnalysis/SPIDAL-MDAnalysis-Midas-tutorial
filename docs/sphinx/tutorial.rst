.. -*- mode: rst; coding: utf-8 -*-

==========
 Tutorial
==========


Preliminary test
================

Provide topology and trajectory files to the psa script as two lists
in a JSON file. Just check that it can process the data ::
  
 $RPDIR/mdanalysis_psa_partial.py --inputfile testcase.json -n 5

This means
- use the test case
- compare the first 5 trajectories against the remaining 5
  trajectories

(The -n (split) argument is important because we are going to use it
to decompose the full distance matrix into sub-matrices. If you just
want to do all-vs-all comparisons, use the ``mdanalysis_psa.py``
script.)

You should see output like ::
   
   Loading paths from JSON file testcase.json
   Processing 10 trajectories.
   Splitting trajectories in two blocks of length 5 and 5
   Calculating D (shape 5 x 5) with 25 entries
   ----------[ TIMING ]--------------------
   load Universes           1.076 s
   PSA distance matrix      2.700 s
   saving output            0.019 s
   ----------------------------------------
   total time               3.795 s
   ----------------------------------------

This indicates that all MDAnalysis parts are working.

Similarly::

   (mdaenv) orbeckst@login2:WORK$ $RPDIR/mdanalysis_psa_partial.py --inputfile testcase.json -n 7
   
   Loading paths from JSON file testcase.json
   Processing 10 trajectories.
   Splitting trajectories in two blocks of length 7 and 3
   Calculating D (shape 7 x 3) with 21 entries
   ----------[ TIMING ]--------------------
   load Universes           0.608 s
   PSA distance matrix      2.325 s
   saving output            0.003 s
   ----------------------------------------
   total time               2.936 s
   ----------------------------------------


Supercomputer environment
=========================

We used `TACC Stampede <https://www.tacc.utexas.edu/stampede/>`_ to
run the calculations on 32 cores but any cluster or even a multicore
workstation might be sufficient.

- Make sure all env vars are set (especially MongoDB,
  :envvar:`RADICAL_PILOT_DBURL`) and password-less ssh works.
- On stampede: Set environment variable  :envvar:`RADICAL_PILOT_PROJECT` to your
  XSEDE allocation::

    export RADICAL_PILOT_PROJECT=TG-xxxxxx

- Activate the *mdaenv* environment.
- You should have the JSON files in the ``WORK`` directory.

Copy the two scripts to the WORK directory (at the moment, this is a
limitation of the scripts to keep them simple) ::

   cd WORK
   cp $RPDIR/{rp_psa.py,mdanalysis_psa_partial.py} .


Radical.pilot script
====================

TODO: describe/show script



Launch pilot jobs
=================


and launch the pilot job::
   
   python rp_psa.py trajectories.json 20 16 spidal_mda_rp_psa 

The ``rp_psa.py`` radical.pilot script takes as input:

- the JSON file with the trajectories (trajectories.json)
- number of trajectories per block (20)
- number of cores to request (16)
- session name (arbitrary string, spidal_mda_rp_psa)


Combine data
============

TODO

Analyze data
============

TODO
