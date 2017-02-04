.. -*- mode: rst; coding: utf-8 -*-

============
 Background
============

We want to perform **Path Similarity Analysis** (PSA) [Seyler2015]_ on
a relatively small data set of 400 molecular dynamics
trajectories. PSA will group similar trajectories with each other and
enables a quantitative comparison. 200 trajectories were produced with
the DIMS MD method [Perilla2011] [Beckstein2009], which is based on
molecular dynamics, and 200 were produced with FRODA [Farrell2010],
which is based on a geometric rigidity decomposition. We want to
understand if these methods produce similar or different transition
pathways, and, perhaps *how* they differ.

Computational cost
==================

PSA requires the calculation of the distances between all paths
:math:`P_i`, the *distance matrix*

.. math::

   D_{ij} = \delta(P_i, P_j)

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

