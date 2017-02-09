.. -*- mode: rst; coding: utf-8 -*-

============
 Background
============

We want to perform **Path Similarity Analysis** (PSA) [Seyler2015]_ on
a relatively small data set of 400 molecular dynamics
trajectories. PSA will group similar trajectories with each other and
enables a quantitative comparison. 200 trajectories were produced with
the DIMS MD method [Perilla2011]_ [Beckstein2009]_, which is based on
molecular dynamics, and 200 were produced with FRODA [Farrell2010]_ ,
which is based on a geometric rigidity decomposition. We want to
understand if these methods produce similar or different transition
pathways, and, perhaps *how* they differ.

The PSA distance metric
=======================

PSA requires the calculation of the distances between all paths
:math:`P_i`, the *distance matrix*

.. math::

   D_{ij} = \delta(P_i, P_j)

The :mod:`MDAnalysis.analysis.psa` module in MDAnalysis_ contains
implementations for various distance functions :math:`\delta(P,
Q)`. Here we are going to use the discrete `Fréchet distance`_ as
implemented in
:func:`MDAnalysis.analysis.psa.discrete_frechet`. Another good choice
is the `Hausdorff distance`_
(:func:`MDAnalysis.analysis.psa.hausdorff`), and it is straightforward
to alter the example code in this tutorial (namely,
:doc:`code/mdanalysis_psa_partial`) to explore the use of a
different distance function.

.. _MDAnalysis: http://mdanalysis.org
.. _Fréchet distance: https://en.wikipedia.org/wiki/Fr%C3%A9chet_distance
.. _Hausdorff distance: https://en.wikipedia.org/wiki/Hausdorff_distance


Computational cost
==================

Because the *path metric* :math:`\delta(P, Q)` is symmetric, the
distance metric is also symmetric. Furthermore, :math:`\delta(P, P) =
0`. Therefore, we need to perform in total

.. math::

   M = \frac{N(N-1)}{2}

calculations. For :math:`N=400`, this means 79,800 individual
trajectory comparisons.


Approach
========

Individual path comparisons are relatively compute intensive (they
scale with :math:`\mathcal{O}(n_i n_j)`, the product of the number of time frames
in each trajectory. But they are all independent and hence a simple
**map-reduce** strategy can be easily implemented whereby individual
jobs compute sub-matrices of the distance matrix (map) and the results
are then combined into the complete solution (reduce).

Block distance matrix
---------------------

We choose a *block size* :math:`w` and decompose the distance matrix
:math:`\mathsf{D} = [D_{ij}],\ 0 \leq i, j < N` into block matrices :math:`\mathsf{D}^{\alpha}`
with :math:`0 \leq \alpha < N/w` and

.. math::

   D^{\alpha}_{u(i),v(j)} = D_{i,j} \quad \text{with}\ u = i - \alpha w,\ v
   = j - \alpha w

(with the appropriate adjustments for blocks at the edges that contain
less than :math:`w` entries).


MDAnalysis
----------

In MDAnalysis each trajectory is loaded into a
:class:`~MDAnalysis.core.universe.Universe`, using a *topology file*
(which contains static information about the atoms, bonds etc) and a
*trajectory file* (which contains the changing coordinates for each
time step)::
  
  import MDAnalysis as mda
  universe = mda.Universe(topology, trajectory)

Because we need to access all trajectory data for our analysis, we can
increase performance by first loading the whole trajectory into memory
[#inmemory]_ with the
:meth:`~MDAnalysis.core.universe.Universe.transfer_to_memory()` method
[#memoryreader]_::

  universe.transfer_to_memory()

We will restrict the calculation of the path distances to a subset of
atoms, namely the "C-alpha" atoms that are distributed along the
protein backbone and report on the larger conformational changes. The
C-alpha atoms can be selected as a
:class:`~MDAnalysis.core.groups.AtomGroup` with the `atom selection
language`_::

  ca = u.select_atoms("name CA")

In order to calculate a pairwise distance between two trajectories we
extract the coordinates of all CA atoms for all time steps into a
numpy array::

  P = u.trajectory.timeseries(asel=ca, format="fac")

With the coordinates for a second trajectory ``Q`` we can then
calculate the discrete Fréchet distance using an recursive dynamic
programmig algorithm [EiterManilla1994]_, implemented in
:func:`MDAnalysis.analysis.psa.discrete_frechet`::

  from MDAnalysis.analysis import psa
  dF = psa.discrete_frechet(P, Q)

``dF`` contains the discrete Fréchet distance :math:`\delta_F(P, Q)`,
a real non-negative number. Because we base the Fréchet distance on
the root mean square distance (RMSD) between the CA coordinates for
two frames in :math:`P` and :math:`Q` as its point-wise metric (see,
for instance, [Seyler2015]_ for more details), :math:`\delta_F(P, Q)`
has the interpretation of the RMSD between the two frames in the two
trajectories that best characterize the difference between the two
trajectories (they form a *Fréchet pair*).

.. Note::

   For macromolecular systems, we typically remove all translational
   and rotational degrees of freedoms for all trajectories by
   superimposing *all* trajectory frames on a single reference
   structure [Seyler2015]_. The superposition can be carried out in a
   pre-processing step using, for instance,
   :class:`MDAnalysis.analysis.align.AlignTraj` or as part of PSA with
   :class:`MDAnalysis.analysis.psa.PSAnalysis`. The trajectories for
   this tutorial were already superimposed appropriately (on the
   "CORE" domain of AdK, as described in more detail in
   [Seyler2015]_.)

Calculating a full Fréchet distance matrix :math:`D_{ij} = \delta(P_i,
P_j)` just requires more book-keeping in order to perform the above
steps for the Cartesian product of all trajectories :math:`P_i \times
P_j`.



Radical.pilot
-------------

We use :mod:`radical.pilot` to generate one *compute unit* for each
block matrix computation. The pilot job distributes the individual
compute units, which includes staging of input trajectories and
retrieval of the output file (the block matrix), as well as running
the MDAnalysis script that performs the calculation of the block
matrix on a compute node.


.. rubric:: Footnotes

.. [#inmemory] Instead of using
   :meth:`~MDAnalysis.core.universe.Universe.transfer_to_memory()` one
   could also simply set the ``in_memory=True`` keyword argument of
   :class:`~MDAnalysis.core.universe.Universe` as shown for `in-memory
   representation of arbitrary trajectories`_. However, here we keep
   the two steps separate for conceptual clarity.

.. [#memoryreader] Loading trajectory data into memory makes use of
   the new :class:`~MDAnalysis.coordinates.memory.MemoryReader`
   functionality in :mod:`MDAnalysis.coordinates.memory`; this will be
   available in the upcoming 0.16.0 release.  The main reason why this
   tutorial is using the current development version of MDAnalysis is
   for using the MemoryReader.


.. _atom selection language:
   http://devdocs.mdanalysis.org/documentation_pages/selections.html

.. _In-memory representation of arbitrary trajectories:
   http:/devdocs.mdanalysis.org/documentation_pages/coordinates/memory.html#in-memory-representation-of-arbitrary-trajectories
