=============================
s4cmb (public version)
=============================

.. image:: https://travis-ci.org/JulienPeloton/s4cmb.svg?branch=master
    :target: https://travis-ci.org/JulienPeloton/s4cmb

.. image:: https://coveralls.io/repos/github/JulienPeloton/s4cmb/badge.svg
    :target: https://coveralls.io/github/JulienPeloton/s4cmb

.. raw:: html

    <img src="https://github.com/JulienPeloton/s4cmb/blob/master/s4cmb/data/intro.png" height="400px">

.. contents:: **Table of Contents**

The package
===============
Systematics For Cosmic Microwave Background (s4cmb), is a package to
study instrumental systematic effects in the context of current and future
Cosmic Microwave Background experiments.

Requirements
===============
The pipeline is mainly written in python and it has the following dependencies:

* numpy, matplotlib
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

Make sure you have correct permissions (otherwise just add --user at the end of the command).
You can also directly use the code by updating manually your PYTHONPATH.
Just add in your bashrc:

::

    s4cmbPATH=/path/to/the/package
    export PYTHONPATH=$PYTHONPATH:$s4cmbPATH

Then run the test suite and the coverage:

::

    ./coverage_and_test.sh

It should print the actual coverage of the test suite, and exit with no errors.

Installation using Docker
===============
Alternatively if you do not want install the package on your computer,
we provide a docker image for s4cmb with always the latest version. Install
docker on your computer, and pull the image:

::

    docker pull julienpeloton/s4cmb:latest

Then create a new container and run an interactive session by just running

::

    docker run -i -t julienpeloton/s4cmb:latest bash

Quick examples
===============
We provide a quick end-to-end example for using the package:

::

    python examples/test/simple_app.py -inifile examples/inifiles/simple_parameters.py -tag test

You can also run it on many processors, using MPI (you will need the package mpi4py):

::

    mpirun -n <nproc> python examples/test/simple_app.py -inifile examples/inifiles/simple_parameters.py -tag test_MPI

where nproc should not be greater than the number of scans to run.
Note that for NERSC users, we also provide a quick submission script for jobs on Cori (see examples/nersc_cori.batch).

s4cmb bootcamp
===============

You can find a bootcamp in two parts (notebooks + examples) at `s4cmb-resources <https://github.com/JulienPeloton/s4cmb-resources>`_.
The goal of this bootcamp is to describe the basic parts of the API, and provide ready-to-use examples (for use on laptop and supercomputer).


TODO
===============

* Add WHWP demodulation module.
* Add correlated noise simulator (and update mapmaking weights).

Main developers
===============
* Julien Peloton (j.peloton at sussex.ac.uk)
* Giulio Fabbian (gfabbian at ias.u-psud.fr)

Support
===============

.. raw:: html

    <img src="https://github.com/JulienPeloton/s4cmb/blob/master/s4cmb/data/LOGO-ERC.jpg" height="200px">
