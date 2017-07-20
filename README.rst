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

Coming soon: dockerfile :-)

Quick examples
===============
You can find notebooks describing how to use basic functionalities of s4cmb
in the folder jupyter_doc.

We also provide a quick end-to-end example for using the package with MPI.
Try to run (you will need the package mpi4py)

::

    mpirun -n <nproc> python examples/simple_app.py -inifile examples/simple_parameters.ini -tag test

where nproc should not be greater than the number of scans to run.
Note that for NERSC users, we also provide a submission script for jobs on Cori (see examples/nersc_cori.batch).

How to build your own s4cmb App?
===============
Let's say we want to build an instrument, a scanning strategy, and scan the sky to obtain
data. Say we also want to inject crosstalk between detectors, and then reconstruct the sky maps with the contamination.

* Step 1 [parameters initialisation]: create a ini file with your parameters. The best is to copy the one provided (examples/simple_parameters.ini) and change the values to yours. Do not forget to update the paths to data!

::

    [s4cmb]
    ## Parameter file for a fake experiment.
    ## Run ID
    tag = gros
    name_instrument = fake

    ...

* Step 2 [start the App]: Create a python script, and import relevant modules

::

    ## python 2/3 compatibility.
    from __future__ import division, absolute_import, print_function

    ## If you want to perform parallel computation.
    from mpi4py import MPI

    ## Import modules and routines from s4cmb.
    import s4cmb

    ...

* Step 3 [tell the App what to read]: link your inifile to your App. For that one we will use the module argparse for example. Also add any useful args you want to pass:

::

    def addargs(parser):
        """ Parse command line arguments for s4cmb """

        ## Defaults args - load instrument, scan and sky parameters
        parser.add_argument(
            '-inifile', dest='inifile',
            required=True,
            help='Configuration file with parameter values.')

        ...

* Step 3 [load background]: Tell the App to load the background (instrument, scan, and so on).

::

    if __name__ == "__main__":
        """
        Launch the pipeline!
        """
        <grab args>

        ## Initialise our input maps.
        sky_in = s4cmb.input_sky.HealpixFitsMap(...)

        ## Initialise our instrument.
        inst = s4cmb.instrument.Hardware(...)

        ## Initialize our scanning strategy and run the scans.
        scan = s4cmb.scanning_strategy.ScanningStrategy(...)
        scan.run()

* Step 4 [perform computations]: Loop over scans, and for each scan do map2tod -> inject crosstalk -> tod2map. Note that the maps are coadded on the fly so that sky_out_tot contains all scans.

::

    for CESnumber in range(scan.nCES):
        tod = s4cmb.tod.TimeOrderedDataPairDiff(...)

        ## Initialise map containers for each processor
        if CESnumber == 0:
            sky_out_tot = s4cmb.tod.OutputSkyMap(...)

        ## Scan input map to get TODs
        d = np.array([
            tod.map2tod(det) for det in range(inst.focal_plane.nbolometer)])

        ## Inject crosstalk
        s4cmb.systematics.inject_crosstalk_inside_SQUID(d, ...)

        ## Project TOD back to maps
        tod.tod2map(d, sky_out_tot)

* Step 5 [write on disk your maps]: We provide some routines to write fits file but feel free to write your routines with your favourite I/O!

::

    s4cmb.xpure.write_maps_a_la_xpure(...)
    s4cmb.xpure.write_weights_a_la_xpure(...)

Et voilà! You can find this complete example in examples/so_crosstalk_app.py.


TODO
===============

* Add the dockerfile.

Main developers
===============
* Julien Peloton (j.peloton at sussex.ac.uk)
* Giulio Fabbian (gfabbian at ias.u-psud.fr)
