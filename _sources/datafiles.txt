.. -*- mode: rst; coding: utf-8 -*-

.. highlight:: bash

============
 Data files
============

In this section you will download all the files that are required for
the tutorial. This includes Python scripts and real-world data files.


.. _install-scripts:

Input files and scripts
-----------------------

Get the latest version of the tutorial::

  git clone https://github.com/MDAnalysis/SPIDAL-MDAnalysis-Midas-tutorial.git

The Python scripts for the tutorial are located in
``SPIDAL-MDAnalysis-Midas-tutorial/rp`` and
``SPIDAL-MDAnalysis-Midas-tutorial/util``. 


Data
----

A real-world data set of about 1.5 GiB size is made available (and
should be downloaded well in advance of this tutorial). These are 400
trajectories of the conformational transition of the enzyme adenylate
kinase [Seyler2015]_, sampled with the dynamic importance sampling
(DIMS) MD method [Beckstein2009]_ [Perilla2011]_ or the FRODA method
[Farrell2010]_.


Download the whole archive as a 1.2 GiB tar.bz2 from
https://becksteinlab.physics.asu.edu/pages/download/SPIDAL-tutorial-data.tar.bz2
and unpack with ::

    tar -jxvf SPIDAL-tutorial-data.tar.bz2

(Alternatively, use the `dropbox PSA folder`_.)

This should create a directory ``tutorial-data``. Do not change the
contents inside the directory because the primitive
:program:`generate_file_list.py` script (see below) has built-in
assumptions on the directory structure.

Put the ``tutorial-data`` directory in the same directory as the
tutorial code directory ``SPIDAL-MDAnalysis-Midas-tutorial``.


Work directory
--------------

Perform all runs in a directory ``WORK`` (here we assume that ``WORK``
in the same directory as the data and the tutorial code)::

   mkdir WORK
   cd WORK


Generate JSON file
------------------

We need a starting JSON file that lists the pairs of (topology,
trajectory) files with absolute paths. This is a short cut for this
tutorial, in general, managing large collections of simulation
trajectories is a difficult topic and tools such as MDSynthesis_ can
be employed.

Generate input file for the RP script :program:`rp/rp_psa.py` with the helper
script :program:`util/generate_file_list.py`::

   cd WORK  # if you are not already in WORK
   RPDIR=$PWD/../SPIDAL-MDAnalysis-Midas-tutorial/rp
   
   $RPDIR/../util/generate_file_list.py -T .. trajectories.json

(Just needs the path to the directory that contains the data and and
output file name.) This creates the list with all 2 x 200 = 400
trajectories.

Also create a smaller list for testing (only 2 x 5 = 10 trajectories)::

   $RPDIR/../util/generate_file_list.py -T .. -e 5 testcase.json

You can have a look at the JSON files with :program:`less` or a text
editor. It is important that the files are listed with absolute paths
to make sure that the radical.pilot script can find them.

.. highlight:: python

.. links

.. _`dropbox PSA folder`:
   https://www.dropbox.com/sh/3sfu6x37lieti26/AAB255qzdUgAQia_XdfoIklCa?dl=0
.. _MDSynthesis: http://mdsynthesis.readthedocs.io/
