.. -*- mode: rst; coding: utf-8 -*-

================
 Available code
================

The tutorial requires three Python scripts, which are listed here for
references and convenience. In order to obtain the scripts, clone the
GitHub repository as described under :ref:`Installing Input Files and
Scripts <install-scripts>`.

Map-reduce scripts
==================

These scripts form the **map-reduce** algorithm. :doc:`rp_psa`
implements the map logic at the radical.pilot level and drives the
computation with :doc:`mdanalysis_psa_partial`. :doc:`psa_reduce`
combines the partial results and performs the analysis.

.. toctree::
   :maxdepth: 1
   :caption: Python scripts:

   rp_psa
   mdanalysis_psa_partial
   psa_reduce


Helper scripts
==============

* :file:`util/generate_file_list.py` is used in this tutorial to
  generate comprehensive input files in JSON format. It is
  specifically written for this tutorial with these specific
  data. Users will need to adapt their own approaches to manage their
  trajectories and communicate topology/trajectory pairs to the
  :doc:`rp_psa` script.
* :file:`rp/mdanalysis_psa.py` is a complete serial implementation of
  the PSA distance matrix calculation; it is provided for comparison
  purposes.
