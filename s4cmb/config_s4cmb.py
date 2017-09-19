#!/usr/bin/python
"""
Module to normalise ini files containing parameter
values for running s4cmb in software mode.
Not required for the API.

Author: Julien Peloton, j.peloton@sussex.ac.uk
"""
import os
import sys
import importlib

def import_string_as_module(fn_full):
    """
    Import module from its name given as a string.

    Parameters
    ----------
    fn_full: string
        Python filename containing parameters that you
        want to import.

    Returns
    ----------
    params: module
        Module containing your parameters.

    Examples
    ----------
    >>> fn_full = 'examples/inifiles/simple_parameters.py'
    >>> params = import_string_as_module(fn_full)
    >>> 'do_pol' in dir(params)
    True

    """
    ## Import parameters from the user parameter file
    fn_short = os.path.basename(fn_full).split('.py')[0]
    sys.path.insert(0, os.path.realpath(
        os.path.join(os.getcwd(), os.path.dirname(fn_full))))
    params = importlib.import_module(fn_short)
    return params


if __name__ == "__main__":
    import doctest
    doctest.testmod()
