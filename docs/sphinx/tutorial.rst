.. -*- mode: rst; coding: utf-8 -*-

==========
 Tutorial
==========

We will first test that MDAnalysis can perform PSA, using the small
test set. We will then look in detail into the radical.pilot script
that implements the map-step. Finally, we conclude with the
reduce-step and analysis of the data. 


Preliminary test
================

Provide topology and trajectory files to the PSA script
:program:`mdanalysis_psa_partial.py` as two lists in a JSON file (see
:ref:`mdanalysis_psa_partial` for more details). Initially we just
want to check that the script can process the data

.. code-block:: bash

   (mdaenv) $ $RPDIR/mdanalysis_psa_partial.py --inputfile testcase.json -n 5

This means

- use the test case
- compare the first 5 trajectories against the remaining 5
  trajectories

The ``-n`` (split) argument is important because we are going to use
it to decompose the full distance matrix into sub-matrices. (If you
just want to do all-vs-all comparisons, use the
:program:`mdanalysis_psa.py` script.)

You should see output like

.. code-block:: none

    Loading paths from JSON file testcase.json
    Processing 10 trajectories.
    Splitting trajectories in two blocks of length 5 and 5
    Calculating D (shape 5 x 5) with 25 entries
    ----------[ TIMING ]--------------------
    load Universes                     0.448 s
    universe.transfer_to_memory()      0.178 s
    PSA distance matrix                1.483 s
    saving output                      0.001 s
    ----------------------------------------
    total time                         2.111 s
    ----------------------------------------

This indicates that all MDAnalysis parts are working.

Similarly:

.. code-block:: bash

   (mdaenv) $ $RPDIR/mdanalysis_psa_partial.py --inputfile testcase.json -n 7

should give something like

.. code-block:: none

    Loading paths from JSON file testcase.json
    Processing 10 trajectories.
    Splitting trajectories in two blocks of length 7 and 3
    Calculating D (shape 7 x 3) with 21 entries
    ----------[ TIMING ]--------------------
    load Universes                     0.422 s
    universe.transfer_to_memory()      0.132 s
    PSA distance matrix                1.436 s
    saving output                      0.001 s
    ----------------------------------------
    total time                         1.991 s
    ----------------------------------------

Supercomputer environment
=========================

We used `TACC Stampede <https://www.tacc.utexas.edu/stampede/>`_ to
run the calculations on 32 cores but any cluster or even a multi core
workstation might be sufficient.

- Make sure all env vars are set (especially MongoDB,
  :envvar:`RADICAL_PILOT_DBURL`) and password-less ssh works.
- On stampede: Set environment variable  :envvar:`RADICAL_PILOT_PROJECT` to your
  XSEDE allocation:

  .. code-block:: bash
  
     export RADICAL_PILOT_PROJECT=TG-xxxxxx

- Activate the *mdaenv* environment.
- You should have the JSON files in the ``WORK`` directory.

Copy the two scripts to the WORK directory (at the moment, this is a
limitation of the scripts to keep them simple)L

.. code-block:: bash

   cd WORK
   cp $RPDIR/{rp_psa.py,mdanalysis_psa_partial.py} .


Map-step: Radical.pilot script
==============================

RP script
---------

The pilot script is :program:`rp/rp_psa.py` (see :ref:`rp_psa` for
details). In the following, important parts of the script are
explained.

First, a session is created and a Pilot Manager is initialized:

.. literalinclude:: /code/rp/rp_psa.py
   :language: python
   :lines: 65-70
   :linenos:
   :lineno-start: 65   

Second, a ComputePilot is described and submitted to the Pilot Manager:

.. literalinclude:: /code/rp/rp_psa.py
   :language: python
   :lines: 73-82
   :linenos: 
   :lineno-start: 73

First, a Unit Manager is created for that session and the Pilot is added:

.. literalinclude:: /code/rp/rp_psa.py
   :language: python
   :lines: 85-88
   :linenos:
   :lineno-start: 85   


The script reads the topology and trajectory files from a JSON file:

.. literalinclude:: /code/rp/rp_psa.py
   :language: python
   :lines: 91-92
   :linenos:
   :lineno-start: 91

The MDAnalysis script is staged

.. literalinclude:: /code/rp/rp_psa.py
   :language: python
   :lines: 99-112
   :linenos:
   :lineno-start: 99

In a loop, a CU is set up for each block matrix --- this is the
**map** step. In particular, the trajectories are partitioned
according to the ``BLOCK_SIZE`` (the parameter :math:`w`) and all
necessary information is written to a JSON file that will be used as
the input for the :program:`mdanalysis_psa_partial.py` script:

.. literalinclude:: /code/rp/rp_psa.py
   :language: python
   :lines: 115-189
   :emphasize-lines: 29-58,70-72
   :linenos:
   :lineno-start: 115   

For the :ref:`reduce step <section-reduce-step>`, all information
about the block matrix (filename and indices in the distance matrix)
are written to a JSON file ("manifest"):

.. literalinclude:: /code/rp/rp_psa.py
   :language: python
   :lines: 193-194
   :linenos: 
   :lineno-start: 193

Finally, the CUs are submitted to execute on the compute resources and
the script waits until they are all complete:

.. literalinclude:: /code/rp/rp_psa.py
   :language: python
   :lines: 197-202
   :linenos:
   :lineno-start: 197



Launch pilot jobs
-----------------

Launch the pilot job from the login node of stampede (or the cluster
where you set up your radical.pilot environment):

.. code-block:: bash
   
   python rp_psa.py trajectories.json 20 16 SPIDAL.001 

The :program:`rp_psa.py` radical.pilot script takes as input (see
:ref:`rp_psa` for details):

- the JSON file with the trajectories (trajectories.json)
- number of trajectories per block (20)
- number of cores to request (16)
- session name (arbitrary string, "SPIDAL.001"); note that you
  need to provide a new session name every time you-rerun the command,
  e.g. "SPIDAL.002", "SPIDAL.003", ...


.. _section-reduce-step:

Reduce-step
===========

Once the pilot jobs has completed and all block matrices have been
computed we combine all data (**reduce**) into the distance matrix
:math:`D_{ij}` and analyze it:

.. code-block:: bash

   psa_reduce.py -t trajectories.json -m manifest.json \
                 -p psa_plot.pdf -o distance_matrix.npy


Combine block matrices
----------------------

The ``manifest.json`` file contains all information that is needed to
re-assemble the distance matrix: for each output file it also stores
the indices of the sub-matrix in the distance matrix and so the full
distance matrix can be built with

.. literalinclude:: /code/rp/psa_reduce.py
   :language: python
   :pyobject: combine
   :emphasize-lines: 4,11,12
   :linenos:

The matrix is written to a numpy file ``distance_matrix.npy`` (and can
be loaded with :func:`numpy.load`).

Analysis
--------

The distance matrix can be clustered to reveal patterns in the
trajectories. Here we use hierarchical clustering with `Ward's linkage
criterion`_ as implemented in :mod:`scipy.cluster.hierarchy`
[Seyler2015]_. The clustering and plotting code was taken from the
:meth:`~MDAnalysis.analysis.psa.PSAnalysis.cluster` and
:meth:`~MDAnalysis.analysis.psa.PSAnalysis.plot` methods of
:class:`~MDAnalysis.analysis.psa.PSAnalysis`. The plotting code is
fairly complicated but the clustering is straightforward:

.. literalinclude:: /code/rp/psa_reduce.py
   :language: python
   :pyobject: cluster
   :linenos: 

In the resulting plot we indicate FRODA trajectories with a heavy
double-bar and DIMS trajectories with a thin single bar.

.. _figure-psa:

.. figure:: /figures/psa_distance_matrix.png
   :alt: Clustered PSA distance matrix

   **PSA of AdK transitions**. Transitions of AdK from the closed to
   the open conformation were generated with the DIMS MD and the FRODA
   method. *Path similarity analysis* (PSA) was performed by
   hierarchically clustering the matrix of Fréchet distances between
   all trajectories, using `Ward's linkage criterion`_.

   +---------+-----------------------+
   | Symbol  | Method                |
   +=========+=======================+
   | `||`    | FRODA                 |
   +---------+-----------------------+
   | `|`     | DIMS                  |
   +---------+-----------------------+

Figure :ref:`PSA of AdK transitions <figure-psa>` shows clearly that DIMS and FRODA transitions
are different from each other. Each method produces transitions that
are characteristic of the method itself. In a comparison of many different
transition path sampling methods it was also found that this pattern
generally holds, which indicates that there is likely no "best" method
yet [Seyler2015]_.



.. _`Ward's linkage criterion`:
   https://en.wikipedia.org/wiki/Ward's_method

   


