v0.x.x
=============
* Current version.
* Add dichroic detector (dichroic branch)
* Add subscans information
* Add beam ellipticity systematics

v0.6.1
=============
* Add Neil small aperture scan
* Add method to perform systematic deprojection (IGQU)
* Add App for differential pointing
* Fix bug in dealing with leap second correction.
* Allow python 3.5 & 3.6 (both code and Travis fixed)
* Allow image in interactive ipython to be shown
* Completely remove weave - never used in practice (use fortran instead for speed-up)

v0.6.0
=============
* Add HWP demodulation (in addition to pair difference) to estimate I, Q, and U from timestreams.
* Allow different nside (input vs output) to be used.

v0.5.3
=============
* Add ACT scanning strategy, and custom scanning strategy.
* Add link to the bootcamp in the description (`s4cmb-resources <https://github.com/JulienPeloton/s4cmb-resources>`_).

v0.5.2
=============
* Improve coverage and coverage tool.
* Simplify the scanning strategy module.
* Add App to study relative gain variation
* Allow to process channel-by-channel in tod2map to save memory!
* Update the default noise level to 5ish uk.arcmin.
* Add routines to simulate gain drifts.

v0.5.1
=============
* Fix bug in bolometer coordinates in the focal plane (it was mixed with pair coordinates).
* Allow routine in HealpixFitsMap to read alms files directly.
* Allow the mapping per detector pair to save memory.
* Remove the notebooks from the repo (now at `s4cmb_notebooks <https://github.com/JulienPeloton/s4cmb_notebooks>`_).
* Add new routines to simulate differential pointing.
* Add new routines to define detectors gain.

v0.5.0
=============
* Include flat sky projection for the output maps.
* Include white noise simulator (time-domain noise).
* Include Dockerfile.

v0.4.0
=============
* Add examples for MPI usage (with examples on Cori, NERSC).
* Add fits interface for dumping maps on disk.
* Remove unused args in several classes and routines.
* Add App for SO, and support for xpure.
* Massive change in the way parsers are handled.
* Add SQUID to SQUID crosstalk.

v0.3.3
=============
* Release the package on pip.

v0.3.1
=============
* Add systematic module, with routine to simulate detector crosstalk.
* Fix bugs in the detector/SQUID labeling (instrument module).

v0.2.0
=============
* Add Fortran scripts and Makefile for the compilation.
* Add detector pointing and TOD modules.
* Doc contains end-to-end example (without systematic effect).
* Fix many bugs related to absolute imports (py 2 vs 3) and failed tests in doctest.

v0.1.0
=============
* First version of the repo (public)
* Add instrument and scanning strategy modules.
* Release public version.
