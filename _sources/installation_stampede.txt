.. -*- mode: rst; coding: utf-8 -*-

===============================================
 Notes on set up for tutorial on TACC Stampede
===============================================

Install miniconda::
  
  wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
  bash Miniconda2-latest-Linux-x86_64.sh

Install in the standard location. Modify bashrc (note: you might want
to undo this later and replace with a more-nuanced approach but for
right now we are trying to get everything working in a simple fashion).

   You may wish to edit your .bashrc or prepend the Miniconda2 install location::

     $ export PATH=/home1/02102/orbeckst/miniconda2/bin:$PATH


Install prereqs::
   
   conda config --add channels MDAnalysis
   conda config --add channels conda-forge
   conda update --yes conda
   conda create --yes -n mdaenv python=2.7 numpy mmtf-python nose=1.3.7 mock sphinx=1.3 six biopython networkx cython joblib griddataformats  
   conda install --yes -n mdaenv matplotlib netcdf4 scikit-learn scipy seaborn coveralls
   conda install --yes -n mdaenv -c biobuilds --yes clustalw=2.1
   
   source activate mdaenv



   git clone --depth=50 https://github.com/MDAnalysis/mdanalysis.git
   cd mdanalysis

   pip install package/
   pip install testsuite/


   pip install radical.pilot

Set up environment.

You need a working MongoDB, for instance, the free plan at
https://mlab.com/ ::
   # export PATH=$HOME/miniconda2/bin:$PATH
   
   source activate mdaenv
   export RADICAL_PILOT_DBURL="mongodb://<user>:<password>@<host>:<port>/<database>"

