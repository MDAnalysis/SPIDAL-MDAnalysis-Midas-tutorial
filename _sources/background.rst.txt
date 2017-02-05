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
to alter the example code from this tutorial to explore the user of a
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

We choose a *block size* :math:`w` and decompose the distance matrix
:math:`\mathsf{D} = [D_{ij}],\ 0 \leq i, j < N` into block matrices.

