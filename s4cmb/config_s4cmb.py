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

def compare_version_number(version, threshold):
    """
    Compare two version numbers.

    Parameters
    ----------
    version: string
        Version of you package x.y.z
    threshold: string
        Threshold version x.y.z

    Returns
    ----------
    result: boolean
        True if your version is higher or equal than the threshold.
        False otherwise.

    Examples
    ----------
    >>> version = '1.10.0'
    >>> threshold = '1.9.1'
    >>> compare_version_number(version, threshold)
    True

    """
    ## If the two versions are equal
    if version == threshold:
        return True

    version_numbers = version.split('.')
    threshold_numbers = threshold.split('.')
    for v, t in zip(version_numbers, threshold_numbers):
        v = int(v)
        t = int(t)
        if v == t:
            continue
        if v > t:
            return True
        if v < t:
            return False
    return True

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
