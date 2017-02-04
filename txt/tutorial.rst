.. -*- mode: rst; coding: utf-8 -*-

=========================================
 Tutorial: MDAnalysis with radical.pilot
=========================================


Setting up packages
===================

conda
-----

Use ``conda`` (with Python 2) with the anaconda_ or miniconda_
 distribution (Python 2.7).

.. _anaconda: https://www.continuum.io/downloads
.. _miniconda: https://conda.io/miniconda.html



MDAnalysis
----------

https://github.com/MDAnalysis/mdanalysis/wiki/Installing%20from%20binary%20packages

Install numpy and other packages::
  
   conda config --add channels MDAnalysis
   conda config --add channels conda-forge
   conda update --yes conda
   conda create --yes -n mdaenv python=2.7 numpy mmtf-python nose=1.3.7 mock sphinx=1.3 six biopython networkx cython joblib griddataformats  
   conda install --yes -n mdaenv matplotlib netcdf4 scikit-learn scipy seaborn coveralls
   conda install --yes -n mdaenv -c biobuilds --yes clustalw=2.1
   
   source activate mdaenv
      
`Install MDAnalysis devel <https://github.com/MDAnalysis/mdanalysis/wiki/Setup-Development-Environment>`_ manually (for the tutorial we will have a snapshot tar ball or perhaps even a conda package)::

  git clone --depth=50 https://github.com/MDAnalysis/mdanalysis.git
  cd mdanalysis
  
Check that you are on develop::
  git branch
Should show::
  * develop

Install development version::
  
  pip install package/
  pip install testsuite/

Check success::

  $ python -c 'import MDAnalysis as mda; print(mda.__version__)'
  0.16.0-dev0
  

Radical Pilot
-------------

http://radicalpilot.readthedocs.io/en/latest/installation.html

In the same env::

  source activate mdaenv
  pip install radical.pilot

Check success::
  $    radicalpilot-version
  0.44
  
You also need a MongoDB. ::

  export RADICAL_PILOT_DBURL="mongodb://user:pass@host:port/dbname"
  export RADICAL_PILOT_DBURL="mongodb://host:port/dbname"

Password-less ssh: http://radicalpilot.readthedocs.io/en/latest/installation.html#setup-ssh-access-to-target-resources

  
.. NOTE:: 

    Many installation problems boil down to one of two causes: an
    Anaconda based Python distribution, or an incompatible version of
    pip/setuptools.

    Many recent systems, specifically in the academic community,
    install Python in its incarnation as Anaconda Distribution. RP is
    not yet able to function in that environment. While support of
    Anaconda is planned in the near future, you will have to revert to
    a ‘normal’ Python distribution to use RP.

Tutorial files
==============

Input files and scripts
-----------------------

Get the latest version of the tutorial::

  git clone https://github.com/Becksteinlab/SPIDAL-MDAnalysis-Midas-tutorial.git



Data
----

Download the whole archive as a 1.2 GiB tar.bz2 from
https://becksteinlab.physics.asu.edu/pages/download/SPIDAL-tutorial-data.tar.bz2
and unpack with ::

    tar -jxvf SPIDAL-tutorial-data.tar.bz2


This should create a directory ``tutorial-data``. Do not change the
contents inside the directory because the ``generate_file_list.py``
script (see below) has built-in assumptions on the structure.

Put the ``tutorial-data`` directory in the same directory as the
tutorial code directory ``SPIDAL-MDAnalysis-Midas-tutorial``.


Work dir
--------

Perform all runs in a directory ``WORK`` (here we assume that ``WORK``
in the same directory as the data and the tutorial code)::

  mkdir WORK
  cd WORK


Generate json file
------------------

We need a starting json file that lists the pairs of (topology,
trajectory) files with absolute paths. This is a short cut for this
tutorial, in general, managing large collections of simulation
trajectories are a difficult topic and tools such as MDSynthesis_ can
be employed.

.. _MDSynthesis: http://mdsynthesis.readthedocs.io/

Generate input file for the RP script ``rp/rp_psa.py`` with the helper
script ``./rp/generate_file_list.py``::

   cd WORK  # if you are not already in WORK
   RPDIR=$PWD/../SPIDAL-MDAnalysis-Midas-tutorial/rp
   
   $RPDIR/generate_file_list.py -T .. trajectories.json

(Just needs the path to the directory that contains the data and and
output file name.) This creates the list with all 2 x 200 = 400
trajectories.

Also create a smaller list for testing (only 2 x 5 = 10 trajectories)::

   $RPDIR/generate_file_list.py -T .. -e 5 testcase.json


Scripts
=======

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


Launch pilot jobs
=================

- Make sure all env vars are set (especially MongoDB,
  :envvar:`RADICAL_PILOT_DBURL`) and password-less ssh works.
- Set environment variable  :envvar:`RADICAL_PILOT_PROJECT` to your
  XSEDE allocation::

    export RADICAL_PILOT_PROJECT=TG-xxxxxx

- Activate the *mdaenv* environment.
- You should have the JSON files in the ``WORK`` directory.

Copy the two scripts to the WORK directory (at the moment, this is a
limitation of the scripts to keep them simple) ::

   cd WORK
   cp $RPDIR/{rp_psa.py,mdanalysis_psa_partial.py} .
   
and launch the pilot job::
   
   python rp_psa.py trajectories.json 20 16 spidal_mda_rp_psa 

The ``rp_psa.py`` radical.pilot script takes as input:

- the JSON file with the trajectories (trajectories.json)
- number of trajectories per block (20)
- number of cores to request (16)
- session name (arbitrary string, spidal_mda_rp_psa)

