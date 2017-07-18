=============================
s4cmb (public version)
=============================

.. image:: https://travis-ci.org/JulienPeloton/s4cmb.svg?branch=master
    :target: https://travis-ci.org/JulienPeloton/s4cmb

.. figure:: https://github.com/JulienPeloton/s4cmb/blob/master/s4cmb/data/intro.png
    :scale: 25 %

The package
===============
Systematics For Cosmic Microwave Background (s4cmb), is a package to
study instrumental systematic effects in the context of current and future
Cosmic Microwave Background experiments.

Requirements
===============
The pipeline is mainly written in python and it has the following dependencies:

* numpy, matplotlib
* h5py (I/O)
* astropy, ephem, pyslalib, healpy (astro libs)
* f2py, weave (interfacing with python)

While we use python 2.7, we try to make it compatible with python 3.x.
If you are using python 3.x and you encounter an error, please open an issue or a
pull request so that we fix it asap.

Some parts of the pipeline are written in C (and compiled on-the-fly via the
package weave), and in Fortran (to come). The latter is interfaced with
python using f2py. The compilation is done usually when you install the
package (see setup.py), but we also provide a Makefile for more
customized compilations (see dir/Makefile).

Installation
===============
You can easily install the package using pip

::

    pip install s4cmb

Otherwise you can fork the repo from the github repository and clone it to your machine.
Use the setup.py for the installation. Just run:

::

    python setup.py install

Make sure you have correct permissions (otherwise just add --user).
You can also directly use the code by updating manually your PYTHONPATH.
Just add in your bashrc:

::

    s4cmbPATH=/path/to/the/package
    export PYTHONPATH=$PYTHONPATH:$s4cmbPATH

Then run the test suite and the coverage:

::

    ./coverage_and_test.sh

It should print the actual coverage of the test suite, and exit with no errors.

Coming soon: dockerfile :-)

Examples
===============
You can find notebooks describing how to use basic functionalities of s4cmb
in the folder jupyter_doc.

We also provide a full example for using the package on clusters.
Try to run (you will need the package mpi4py)

::

    mpirun -n <nproc> python examples/simple_app.py -inifile examples/simple_parameters.ini

where nproc should not be greater than the number of scans to run.
Note that for NERSC users, we also provide a submission script for jobs on Cori
 (see examples/nersc_cori.batch).

TODO
===============

* Add the dockerfile.

Main developers
===============
* Julien Peloton (j.peloton at sussex.ac.uk)
* Giulio Fabbian (gfabbian at ias.u-psud.fr)
