=============================
s4cmb
=============================

.. image:: https://github.com/JulienPeloton/s4cmb/workflows/s4cmb/badge.svg
    :target: https://github.com/JulienPeloton/s4cmb/actions?query=workflow%3As4cmb

.. image:: https://github.com/JulienPeloton/s4cmb/workflows/PEP8/badge.svg
    :target: https://github.com/JulienPeloton/s4cmb/actions?query=workflow%3APEP8

.. image:: https://coveralls.io/repos/github/JulienPeloton/s4cmb/badge.svg?branch=master
    :target: https://coveralls.io/github/JulienPeloton/s4cmb?branch=master

.. raw:: html

    <img src="https://github.com/JulienPeloton/s4cmb/blob/master/s4cmb/data/intro.png" height="400px">

.. contents:: **Table of Contents**

The package
===============
Systematics For Cosmic Microwave Background (s4cmb), is a package to
study instrumental systematic effects in the context of current and future
Cosmic Microwave Background experiments. Currently accessible:

* Electrical crosstalk in the multiplexed readout.
* Relative gain-calibration uncertainty between the two detectors in a focal plane pixel.
* Time drift of the gains between two consecutive calibration measurements.
* Differential pointing between the two detectors in a pixel.
* ... more to come!

Requirements
===============
The pipeline is mainly written in python (>= 3.6) and it has the following dependencies (see requirements.txt):

* numpy, matplotlib
* astropy, ephem, pyslalib, healpy (astro libs)
* f2py (interfacing with python)

Some parts of the pipeline are written in Fortran which is interfaced with
python using f2py. The compilation is done usually when you install the
package (see setup.py), but we also provide a Makefile for more
customized compilations (see the Makefile in s4cmb).

Installation
===============

The type of installation depends on what you want to do with the code:
just using it or also developing it?

**I just want to use the code:**

You can easily install the package using pip

::

    pip install s4cmb

**In addition to use the code, I want to be a developer:**

The best is to fork the repo from this github repository to your account and clone it to your machine.
Once you have the repo cloned on your machine, use the makefile to compile the source

::

    cd /path/to/s4cmb
    pip install -r requirements.txt
    make

Do not forget to update your PYTHONPATH. Just add in your bashrc:

::

    s4cmbPATH=/path/to/the/s4cmb
    export PYTHONPATH=$PYTHONPATH:$s4cmbPATH

Then run the test suite and the coverage:

::

    ./coverage_and_test.sh

It should print the actual coverage of the test suite, and exit with no errors.

Installation and usage at NERSC
===============

Again, you can easily install the package using pip

::

    pip install s4cmb --user

Alternatively, if you want to do dev at NERSC and do a manual installation, it's better to keep most of your packages under Anaconda.
I recommend to have a look first at the `NERSC page <https://www.nersc.gov/users/data-analytics/data-analytics-2/python/anaconda-python/>`_ describing how to use it.

The installation of s4cmb can be done in few steps:

* Clone the repo somewhere in your $HOME
* Install dependencies (see requirements.txt) using Anaconda
* Compile the source (using make in /path/s4cmb)

Working with Docker
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
* Julien Peloton (peloton at lal.in2p3.fr)
* Giulio Fabbian (g.fabbian at sussex.ac.uk)

Thanks to
===============
* @ngoecknerwald: original author for a large part of the scanning strategy module.
* @giuspugl, @dpole, @joydidier, and all `contributors <https://github.com/JulienPeloton/s4cmb/graphs/contributors>`_ for all valuable comments, tests, and feedbacks!

In the literature
===============

s4cmb has been used in several research projects:

* Instrumental systematics biases in CMB lensing reconstruction: a simulation-based assessment (`2011.13910 <https://arxiv.org/abs/2011.13910>`_)
* Development of Calibration Strategies for the Simons Observatory (`1810.04633 <https://arxiv.org/abs/1810.04633>`_)
* Studies of Systematic Uncertainties for Simons Observatory: Detector Array Effects (`1808.10491 <https://arxiv.org/abs/1808.10491>`_)
* Studies of Systematic Uncertainties for Simons Observatory: Polarization Modulator Related Effects (`1808.07442 <https://arxiv.org/abs/1808.07442>`_)
* Iterative map-making with two-level preconditioning for polarized Cosmic Microwave Background data sets (`1801.08937 <https://arxiv.org/abs/1801.08937>`_)

Support
===============

.. raw:: html

    <img src="https://github.com/JulienPeloton/s4cmb/blob/master/s4cmb/data/LOGO-ERC.jpg" height="200px">
