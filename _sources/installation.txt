.. -*- mode: rst; coding: utf-8 -*-
.. highlight:: bash

==============
 Installation
==============

The tutorial will demonstrate how one can use the Midas framework,
namely, the `Pilot job`_ abstractions in `radical.pilot`_ in order to
easily parallelize a large number of calculations with the MDAnalysis_
library. First we need to install **MDAnalysis** and
**radical.pilot**.

We will use the conda_ package manager and set up a virtual
environment to keep the computational environment (named *mdaenv*)
clean for the tutorial.

.. _conda: https://conda.io/docs/intro.html
.. _Pilot job: https://en.wikipedia.org/wiki/Pilot_job
.. _radical.pilot: https://radicalpilot.readthedocs.io/
.. _MDAnalysis: http://mdanalysis.org

.. note:: You need to set up working environments on any machine where
          you want to *launch* jobs (e.g., a cluster head node, a
          workstation, or your laptop) and *run* jobs (i.e., typically
          the compute nodes of your cluster).


Installing MDAnalysis
=====================

MDAnalysis is, in principle, full installable with `conda`_, but in
this tutorial we will use the latest development version and
therefore, we will first install pre-requisites with ``conda`` and
then install MDAnalysis from source.

conda
-----

Install the anaconda_ or miniconda_ distribution (Python 2.7) by
following the links and the instructions. Install Python 2.7.

.. _anaconda: https://www.continuum.io/downloads
.. _miniconda: https://conda.io/miniconda.html

For example, on Linux::

  wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
  bash Miniconda2-latest-Linux-x86_64.sh

Accept the modification of your start-up files; if you later do not
want to use the anaconda distribution, remove a line like ::

  export PATH=/home1/02102/orbeckst/miniconda2/bin:$PATH

from your :file:`~/.bashrc` file.


Dependencies of MDAnalysis
--------------------------

Install numpy and other packages for a feature-complete ("full")
installation::
  
   conda config --add channels MDAnalysis
   conda config --add channels conda-forge
   conda update --yes conda
   conda create --yes -n mdaenv python=2.7 numpy mmtf-python \
              nose=1.3.7 mock sphinx=1.3 six biopython networkx \
	      cython joblib griddataformats  
   conda install --yes -n mdaenv matplotlib netcdf4 scikit-learn \
              scipy seaborn coveralls
   conda install --yes -n mdaenv -c biobuilds --yes clustalw=2.1

Activate the installation (has to be done in every shell)::

   source activate mdaenv

MDAnalysis source installation
------------------------------

`Install MDAnalysis devel
<https://github.com/MDAnalysis/mdanalysis/wiki/Setup-Development-Environment>`_
manually::

  git clone --depth=50 https://github.com/MDAnalysis/mdanalysis.git
  cd mdanalysis
  
Check that you are on the development branch::

  git branch

Should show::

  * develop

Install development version::
  
  pip install package/
  pip install testsuite/

Check success::

  $ python -c 'import MDAnalysis as mda; print(mda.__version__)'
  0.16.0-dev0
  

Installing Radical Pilot
========================

We follow the `radical.pilot installation instructions`_

.. _`radical.pilot installation instructions`:
   http://radicalpilot.readthedocs.io/en/latest/installation.html

Install in the same virtual environment [#venv]_::

  pip install radical.pilot

Check success::

  $    radicalpilot-version
  0.44
  
You also `need a MongoDB`_ instance and set the environment variable
:envvar:`RADICAL_PILOT_DBURL` either with full username and password
information ::

  export RADICAL_PILOT_DBURL="mongodb://user:pass@host:port/dbname"

or, for open installations, just ::

  export RADICAL_PILOT_DBURL="mongodb://host:port/dbname"

Note that you can set up a free MongoDB at http://mlab.com.

.. _need a MongoDB:
   https://radicalpilot.readthedocs.io/en/latest/installation.html#mongodb-service

Finally, set up `password-less ssh`_ (or password-less gsissh with
certificates).

.. _password-less ssh:
  http://radicalpilot.readthedocs.io/en/latest/installation.html#setup-ssh-access-to-target-resources

 
Setting up the environment
==========================

Although you don't need to install the software every time, you do
need to ensure that you always work in the configured virtual
environment *mdaenv* [#venv]_ and that you set set the environment
variable for the MongoDB::

  source activate mdaenv
  export RADICAL_PILOT_DBURL="mongodb://user:pass@host:port/dbname"  


.. rubric:: Footnotes

.. [#venv] Do not forget to activate the *mdaenv* environment whenever
           you open a new terminal::
	  
            source activate mdaenv

           Otherwise, you will probably find that scripts cannot find
           MDAnalysis or radical.pilot.


.. highlight:: python
