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

MDAnalysis is full installable with `conda`_ (but if you want to, you
can also install it directly from source [#devinstall]_).

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

  export PATH=/home/USERNAME/miniconda2/bin:$PATH

from your :file:`~/.bashrc` file.


MDAnalysis and dependencies
---------------------------

The current MDAnalysis 0.16.x only fully supports Python 2.7 so we
will create a Python 2.7 environment. By installing the *MDAnalysis*
and *MDAnalysisTests* packages we automatically install all
dependencies as well (such as numpy and scipy as well as various
plotting and data processing libraries that will be used for plotting
the results) [#clustalw]_::
  
   conda config --add channels MDAnalysis
   conda config --add channels conda-forge
   conda update --yes conda
   
   conda create --yes -n mdaenv python=2.7
   conda install --yes -n mdaenv mdanalysis mdanalysistests

Activate the installation (has to be done in every shell)::

   source activate mdaenv

Check success [#prompt]_::

  (mdaenv) $ python -c 'import MDAnalysis as mda; print(mda.__version__)'
  0.16.0
  

Installing Radical Pilot
========================

We follow the `radical.pilot installation instructions`_ — see there
for details. The following are all the required steps but they leave
out background details or trouble shooting tips. If something is not
working, please see `radical.pilot installation instructions`_ or, if
you are at a workshop, ask an instructor.

Install in the same virtual environment [#venv]_ as MDAnalysis
[#prompt]_::

  (mdaenv) $ pip install radical.pilot

Check success [#prompt]_::

  (mdaenv) $ radicalpilot-version
  0.45.1
  
You also `need a MongoDB`_ instance and set the environment variable
:envvar:`RADICAL_PILOT_DBURL` either with full username and password
information ::

  export RADICAL_PILOT_DBURL="mongodb://user:pass@host:port/dbname"

or, for open installations, just ::

  export RADICAL_PILOT_DBURL="mongodb://host:port/dbname"

Note that you can set up a free MongoDB at http://mlab.com.

If you are doing this tutorial as part of a workshop, the organizers
will provide you with a working MongoDB installation and instruct you
what to put into :envvar:`RADICAL_PILOT_DBURL`.

Finally, set up `password-less ssh`_ (or password-less gsissh with
certificates).

.. _`radical.pilot installation instructions`:
   http://radicalpilot.readthedocs.io/en/latest/installation.html

.. _need a MongoDB:
   https://radicalpilot.readthedocs.io/en/latest/installation.html#mongodb-service

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

.. [#devinstall] If for any reason you want to perform a **MDAnalysis
   source installation**, follow these steps to `install MDAnalysis
   devel
   <https://github.com/MDAnalysis/mdanalysis/wiki/Setup-Development-Environment>`_
   manually.

   The following assumes that you are working in a virtual environment
   or that you do a user installation (``pip install --user`` but in
   the following the ``--user`` flag is omitted). It also assumes that
   you have a complete compiler tool chain installed.

   First get pre-requisites::

     pip install numpy cython

   Get the latest source code from the source code repository (you
   need to have git_ installed)::

     git clone --depth=50 https://github.com/MDAnalysis/mdanalysis.git
     cd mdanalysis

   Check that you are on the development branch::

     git branch

   Should show::

     * develop

   Install development version::

     pip install package/
     pip install testsuite/

   Check as above that you can ``import MDAnalysis`` (see main text)::

     $ python -c 'import MDAnalysis as mda; print(mda.__version__)'
     0.16.1-dev

   This should show a release number with the *-dev* suffix to
   indicate a development version.

   .. _git: https://git-scm.com/

.. [#clustalw] Optionally, for a full feature environment you may also
   install :program:`clustalw2` (for sequence alignments with
   :func:`MDAnalysis.analysis.align.fasta2select`)::

      conda install --yes -n mdaenv -c biobuilds --yes clustalw=2.1

   This is not required for this tutorial but if you start using
   MDAnalysis for other projects, you might want to be able to use all
   features without having to think about installing additional
   optional dependencies later.

.. [#prompt] In the following, the shell's prompt is shown as
             ``(mdaenv) $`` and should *not* be typed. It is supposed
             to remind you that you *must be in the virtual
             environment* [#venv]_. Only type what follows after the
             prompt. If the commands give any output, it is shown on
             the lines following the input.

.. [#venv] Do not forget to activate the *mdaenv* environment whenever
           you open a new terminal::
	  
            source activate mdaenv

           Otherwise, you will probably find that scripts cannot find
           MDAnalysis or radical.pilot.


.. highlight:: python



