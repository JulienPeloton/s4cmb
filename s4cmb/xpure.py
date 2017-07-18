#!/usr/bin/python
"""
Script to generate normalized inputs for the software x2pure
(pure pseudo-spectrum estimator).
For > 99.9 percent of the population, this module is useless and you do not
have to use it.

Author: Julien Peloton, j.peloton@sussex.ac.uk
"""
from __future__ import division, absolute_import, print_function

import os
import numpy as np

from s4cmb.input_sky import write_healpix_cmbmap

def safe_mkdir(path, verbose=False):
    """
    Create a path and catch the race condition between path exists and mkdir.

    Parameters
    ----------
    path : string
        Name of the folder to be created.
    verbose : bool
        If True, print messages about the status.

    Examples
    ----------
    Folders are created
    >>> safe_mkdir('toto/titi/tutu', verbose=True)

    Folders aren't created because they already exist
    >>> safe_mkdir('toto/titi/tutu', verbose=True)
    Folders toto/titi/tutu already exist. Not created.

    >>> os.removedirs('toto/titi/tutu')
    """
    abspath = os.path.abspath(path)
    if not os.path.exists(abspath):
        try:
            os.makedirs(abspath)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
    else:
        if verbose:
            print("Folders {} already exist. Not created.".format(path))

def qu_weight_mineig(cc, cs, ss, epsilon=0., verbose=False):
    """
    Create a weight map using the smallest eigenvalue
    of the polarization matrix.

    Parameters
    ----------
    cc : 1d array
        Projected noise weighted cosine**2 per-pixel.
    ss : 1d array
        Projected noise weighted sine**2 per-pixel.
    cs : 1d array
        Projected noise weighted cosine * sine per-pixel.
    epsilon : float
        Threshold for selecting the pixels. 0 <= epsilon < 1/4.
        The higher the more selective.

    Returns
    ----------
    weights : 1d array
        Polarisation weights per-pixel.

    """
    trace = cc + ss
    trace2 = trace * trace
    det = (cc * ss - cs * cs)

    val2 = trace2 - 4 * det
    valid = (val2 > 0.0)
    val = np.zeros_like(val2)
    val[valid] = np.sqrt(val2[valid])

    weight = np.zeros_like(trace)
    lambda_minus = (trace - val) / 2

    if verbose:
        print('criterion is', epsilon, '< det < 1/4 (epsilon= 0. by default)')

    valid = (lambda_minus > (trace - np.sqrt(trace2 - 4*epsilon * trace2)) / 2)

    if verbose:
        valid3 = [x for x in valid if x is True]
        print('number of pixels kept:', len(valid3), '/', np.sum(trace > 0))
        print('Percentage cut: {:3f} %%'.format(
            (1. - float(len(valid3)) / np.sum(trace > 0)) * 100.))

    weight[valid] = lambda_minus[valid]
    return weight

def write_maps_a_la_xpure(OutputSkyMap, name_out, output_path):
    """
    Grab sky data from OutputSkyMap and write them into files readable by
    the software xpure.

    Parameters
    ----------
    OutputSkyMap : OutputSkyMap instance
        Instance of OutputSkyMap containing map data.
    name_out : string
        File name (.fits) where to store the data.
    output_path : string
        Folder where to put the data.

    """
    safe_mkdir(os.path.join(output_path, name_out))

    fits_I = np.zeros((12 * OutputSkyMap.nside * OutputSkyMap.nside))
    fits_Q = np.zeros((12 * OutputSkyMap.nside * OutputSkyMap.nside))
    fits_U = np.zeros((12 * OutputSkyMap.nside * OutputSkyMap.nside))

    I, Q, U = OutputSkyMap.get_IQU()
    fits_I[OutputSkyMap.obspix] = I
    fits_Q[OutputSkyMap.obspix] = Q
    fits_U[OutputSkyMap.obspix] = U

    full_path = os.path.join(
        output_path, name_out, 'IQU_{}.fits'.format(name_out))

    write_healpix_cmbmap(full_path,
                         data=[fits_I, fits_Q, fits_U],
                         nside=OutputSkyMap.nside,
                         fits_IDL=False,
                         coord=None,
                         colnames=['I', 'Q', 'U'],
                         partial=False,
                         nest=False)

    del fits_I, fits_Q, fits_U

def write_weights_a_la_xpure(OutputSkyMap, name_out, output_path, epsilon,
                             HWP=False):
    """
    Grab weight and mask from OutputSkyMap and write them into files
    readable by the software xpure.

    Parameters
    ----------
    OutputSkyMap : OutputSkyMap instance
        Instance of OutputSkyMap containing map data.
    name_out : string
        File name (.fits) where to store the data.
    output_path : string
        Folder where to put the data.
    epsilon : float
        Threshold for selecting the pixels. 0 <= epsilon < 1/4.
        The higher the more selective.
    HWP : bool
        If True, use demodulation syntax for the weights (w0, w4).
        Default is False (pair difference syntax: w, cc, ss, cs)

    """
    safe_mkdir(os.path.join(output_path, name_out))

    ## Intensity
    weight = np.zeros((12 * OutputSkyMap.nside * OutputSkyMap.nside))

    if not HWP:
        weight[OutputSkyMap.obspix] = OutputSkyMap.w
    else:
        weight[OutputSkyMap.obspix] = OutputSkyMap.w0

    full_path = os.path.join(
        output_path, name_out, 'Iw_{}.fits'.format(name_out))
    write_healpix_cmbmap(full_path,
                         data=weight,
                         nside=OutputSkyMap.nside,
                         fits_IDL=False,
                         coord=None,
                         colnames=['I'],
                         partial=False,
                         nest=False)

    binary = np.zeros((12 * OutputSkyMap.nside * OutputSkyMap.nside))
    mask = np.where((weight > 0))[0]
    binary[mask] = 1.0
    full_path = os.path.join(
        output_path, name_out, 'Iw_{}_norm.fits'.format(name_out))
    write_healpix_cmbmap(full_path,
                         data=binary,
                         nside=OutputSkyMap.nside,
                         fits_IDL=False,
                         coord=None,
                         colnames=['I'],
                         partial=False,
                         nest=False)

    ## Polarisation
    weight = np.zeros((12 * OutputSkyMap.nside * OutputSkyMap.nside))
    if not HWP:
        weight[OutputSkyMap.obspix] = qu_weight_mineig(
            OutputSkyMap.cc, OutputSkyMap.cs, OutputSkyMap.ss, epsilon)
    else:
        weight[OutputSkyMap.obspix] = OutputSkyMap.w4

    full_path = os.path.join(
        output_path, name_out, 'Pw_{}.fits'.format(name_out))
    write_healpix_cmbmap(full_path,
                         data=weight,
                         nside=OutputSkyMap.nside,
                         fits_IDL=False,
                         coord=None,
                         colnames=['P'],
                         partial=False,
                         nest=False)

    binary = np.zeros((12 * OutputSkyMap.nside * OutputSkyMap.nside))
    mask = np.where((weight > 0))[0]
    binary[mask] = 1.0
    full_path = os.path.join(
        output_path, name_out, 'Pw_{}_norm.fits'.format(name_out))
    write_healpix_cmbmap(full_path,
                         data=binary,
                         nside=OutputSkyMap.nside,
                         fits_IDL=False,
                         coord=None,
                         colnames=['P'],
                         partial=False,
                         nest=False)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
