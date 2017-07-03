s4cmb (public version)
==

### The package
Systematics For Cosmic Microwave Background (s4cmb), is a package to
study instrumental systematic effects in the context of current and future
Cosmic Microwave Background experiments.

### Requirements
The pipeline is mainly written in python and it has the following dependencies:
* numpy, scipy, matplotlib, etc.
* f2py (for fortran to python interfacing)
* mpi4py (for parallel computing)
* ...

While we use python 2.7, we try to make it compatible with python 3.x.
If you are using python 3.x and you encounter an error, please open an issue or a
pull request so that we fix it asap.

Some parts of the pipeline are written in C (and compiled on-the-fly via the
package weave), and in Fortran. The latter is interfaced with python using f2py.
The compilation is done usually when you install the package (see setup.py), but
we also provide a Makefile for more customized compilations (see <dir>/Makefile).

### Installation
You can easily install the package using pip
```bash
pip install s4cmb
```

Otherwise, you can clone the repo from the github repository and
use the setup.py for the installation. Just run:
```bash
python setup.py install
```
Make sure you have correct permissions (otherwise just add --user).
You can also directly use the code by updating manually your PYTHONPATH.
Just add in your bashrc:
```bash
s4cmbPATH=/path/to/the/package
export PYTHONPATH=$PYTHONPATH:$s4cmbPATH
```

Alternatively if you do not want install the package on your computer,
you can also use the dockerfile provided to create a docker image:
```bash
docker build -t s4cmb .
docker run -it s4cmb
```

### Main developers
* Julien Peloton (j.peloton at sussex.ac.uk)
* Giulio Fabbian (gfabbian at ias.u-psud.fr)
