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

Coming soon: dockerfile :-)

Examples
===============
You can find notebooks describing how to use basic functionalities of s4cmb
in the folder jupyter_doc.

TODO
===============

* Add the dockerfile.

Main developers
===============
* Julien Peloton (j.peloton at sussex.ac.uk)
* Giulio Fabbian (gfabbian at ias.u-psud.fr)
