.. -*- mode: rst; coding: utf-8 -*-

.. _rp_psa:

===========
 rp_psa.py
===========

The :file:`rp/rp_psa.py` script sets up and executes the radical.pilot
job. It implements the *map* logic and splits the distance matrix
computation into sub-blocks (one block per compute unit). It uses
:doc:`mdanalysis_psa_partial` to perform the calculation. Importantly,
it generates a JSON file with the "key/value" pairs (indices into the
full PSA distance matrix/filenames of the sub-block matrices) that are
necessary for the *reduce* step that is performed with
:doc:`psa_reduce`.

.. literalinclude:: rp/rp_psa.py
   :language: python
   :linenos:
