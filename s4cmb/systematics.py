#!/usr/bin/python
# Copyright (c) 2016-2021 Julien Peloton, Giulio Fabbian.
#
# This file is part of s4cmb
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
"""
Module to handle instrument systematics.

Author: Julien Peloton, peloton@lal.in2p3.fr
        Giulio Fabbian, g.fabbian@sussex.ac.uk
"""
from __future__ import division, absolute_import, print_function

import numpy as np
from scipy import signal

from s4cmb.systematics_f import systematics_f
from s4cmb.instrument import construct_beammap, coordinates_on_grid, gauss2d

arcsecond2rad = np.pi / 180.0 / 3600.0
arcmin2rad = np.pi / 180.0 / 60.0
deg2rad = np.pi / 180.0


def inject_crosstalk_inside_SQUID(
    bolo_data,
    squid_ids,
    bolo_ids,
    mu=-3.0,
    sigma=1.0,
    radius=1,
    beta=2,
    seed=5438765,
    new_array=None,
    language="python",
):
    """
    Introduce leakage between neighboring bolometers within a SQUID.
    You have to provide the list of bolometers, to which SQUID they
    belong, and the index of bolometers within each SQUID.

    Parameters
    ----------
    bolo_data : list of arrays
        Contains the time-ordered data of bolometers.
    squid_ids : list of strings
        Contains the SQUID id for each bolometer.
        Should have the same shape[0] as bolo_data.
    bolo_ids : list of strings
        Contains the positions of bolometers within their SQUID.
        Should have the same shape[0] as bolo_data.
    mu : float, optional
        Mean of the Gaussian used to generate the level of leakage,
        in percent. E.g. mu=1.0 means the leakage coefficients will be
        centred around 1%.
    sigma : float, optional
        Width of the Gaussian used to generate the level of leakage,
        in percent. E.g. sigma=1.0 means the leakage coefficients will be
        of order X \pm 1%.
    radius : int, optional
        Controls the number of bolometers talking within a SQUID.
        radius=n means bolometers [N-n, N-n+1, ..., N, ..., N+n-1, N+n]
        will be crosstalking. `radius` must be at most equal to the number of
        bolometers within a SQUID.
    beta : int, optional
        Exponent controling the attenuation. The idea is that crosstalk depends
        on the separation length. Bolometer $b_i$ crosstalking with bolometers
        $b_{i \pm d}$ ($d \leq FF$) with leakage coefficients
        $\alpha_{i \pm d}$ will see his timestream modified as:
            $b_i = b_i + \dfrac{\alpha_{i \pm d}}{d^\beta}  b_{i \pm d}$.
        Assuming some fixed frequency spacing, hardware consideration
        suggests $\beta = 2$.
    seed : int, optional
        Control the random seed used to generate leakage coefficients.
    new_array : None or bolo_data-like array, optional
        If not None, return a new array of timestreams with the modifications.
        Modify bolo_data directly otherwise. Default is None.
    language : string, optional
        Language to perform computations to be chosen in ['python', 'fortran'].
        Default is python. [fortran is very slow...]

    Example
    ----------
    >>> from s4cmb.tod import load_fake_instrument
    >>> from s4cmb.tod import TimeOrderedDataPairDiff
    >>> inst, scan, sky_in = load_fake_instrument()
    >>> tod = TimeOrderedDataPairDiff(inst, scan, sky_in, CESnumber=0)
    >>> d = np.array([tod.map2tod(det) for det in range(2 * tod.npair)])
    >>> print(round(d[0][0], 3))
    40.95

    Inject crosstalk between neighbour bolometers (radius=1)
    >>> squid_ids = inst.focal_plane.get_indices('Sq')
    >>> bolo_ids = inst.focal_plane.bolo_index_in_squid
    >>> inject_crosstalk_inside_SQUID(d, squid_ids, bolo_ids, radius=1)
    >>> print(round(d[0][0], 3))
    39.561

    One can also keep the original timestreams, and return a new array
    containing modified timestreams.
    >>> d = np.array([tod.map2tod(det) for det in range(2 * tod.npair)])
    >>> d_new = np.zeros_like(d)
    >>> inject_crosstalk_inside_SQUID(d, squid_ids, bolo_ids, radius=1,
    ...     new_array=d_new, language='python')
    >>> print(round(d[0][0], 3), round(d_new[0][0], 3))
    ... #doctest: +NORMALIZE_WHITESPACE
    40.95 39.561

    For large number of bolometers per SQUID, you would prefer fortran
    to python to perform the loops. Choose python otherwise.
    """
    # Make mu and sigma unitless (user provides percentage)
    mu = mu / 100.0
    sigma = sigma / 100.0

    combs = {}
    for bolo in range(len(bolo_data)):
        sq = squid_ids[bolo]
        if sq not in combs:
            combs[sq] = []
        combs[sq].append((bolo_ids[bolo], bolo))

    tsout = 0.0 + bolo_data
    # How much to leak from one bolometer to its neighboring channels
    state = np.random.RandomState(seed)
    cross_amp = state.normal(mu, sigma, len(bolo_data))

    if language == "python":
        for sq in combs:
            for ch, i in combs[sq]:
                for ch2, i2 in combs[sq]:
                    separation_length = abs(ch - ch2)
                    if separation_length > 0 and separation_length <= radius:
                        tsout[i] += (
                            cross_amp[i2] / separation_length ** beta * tsout[i2]
                        )

    elif language == "fortran":
        # F2PY convention
        tsout = np.array(tsout, order="F")
        for sq in combs:
            local_indices = np.array(combs[sq]).flatten()[::2]
            global_indices = np.array(combs[sq]).flatten()[1::2]
            systematics_f.inject_crosstalk_inside_squid_f(
                tsout,
                local_indices,
                global_indices,
                radius,
                cross_amp,
                beta,
                len(local_indices),
                len(bolo_data),
                len(bolo_data[0]),
            )

    if new_array is not None:
        new_array[:] = tsout
    else:
        bolo_data[:] = tsout


def inject_crosstalk_SQUID_to_SQUID(
    bolo_data,
    squid_ids,
    bolo_ids,
    mu=-3.0,
    sigma=1.0,
    squid_attenuation=100.0,
    seed=5438765,
    new_array=None,
):
    """
    Introduce leakage between bolometers from different SQUIDs.
    You have to provide the list of bolometers, to which SQUID they
    belong, and the index of bolometers within each SQUID.

    Parameters
    ----------
    bolo_data : list of arrays
        Contains the time-ordered data of bolometers.
    squid_ids : list of strings
        Contains the SQUID id for each bolometer.
        Should have the same shape[0] as bolo_data.
    bolo_ids : list of strings
        Contains the positions of bolometers within their SQUID.
        Should have the same shape[0] as bolo_data.
    mu : float, optional
        Mean of the Gaussian used to generate the level of leakage,
        in percent. E.g. mu=1.0 means the leakage coefficients will be
        of order 1%.
    sigma : float, optional
        Width of the Gaussian used to generate the level of leakage,
        in percent. E.g. sigma=1.0 means the leakage coefficients will be
        of order X \pm 1%.
    squid_attenuation : int, optional
        Value controling the attenuation due to distance between SQUIDS.
        The idea is that crosstalk depends on the separation
        length between SQUIDs. Bolometer $b_i$ crosstalking with bolometers
        $b_{j}$ with leakage coefficients
        $\alpha_{j}$ will see his timestream modified as:
            $b_i = b_i + \dfrac{\alpha_{j}}{squid_attenuation}  b_{j}$.
        Default is 100.
    seed : int, optional
        Control the random seed used to generate leakage coefficients.
    new_array : None or bolo_data-like array, optional
        If not None, return a new array of timestreams with the modifications.
        Modify bolo_data directly otherwise. Default is None.

    Example
    ----------
    >>> from s4cmb.tod import load_fake_instrument
    >>> from s4cmb.tod import TimeOrderedDataPairDiff
    >>> inst, scan, sky_in = load_fake_instrument(nsquid_per_mux=4)
    >>> tod = TimeOrderedDataPairDiff(inst, scan, sky_in, CESnumber=0)
    >>> d = np.array([tod.map2tod(det) for det in range(2 * tod.npair)])
    >>> print(round(d[0][0], 3))
    40.948

    Inject crosstalk between bolometers in different SQUIDs
    >>> squid_ids = inst.focal_plane.get_indices('Sq')
    >>> bolo_ids = inst.focal_plane.bolo_index_in_squid
    >>> inject_crosstalk_SQUID_to_SQUID(d, squid_ids, bolo_ids,
    ...     squid_attenuation=100.)
    >>> print(round(d[0][0], 3))
    40.716

    One can also keep the original timestreams, and return a new array
    containing modified timestreams.
    >>> d = np.array([tod.map2tod(det) for det in range(2 * tod.npair)])
    >>> d_new = np.zeros_like(d)
    >>> inject_crosstalk_SQUID_to_SQUID(d, squid_ids, bolo_ids,
    ...     squid_attenuation=100., new_array=d_new)
    >>> print(round(d[0][0], 3), round(d_new[0][0], 3))
    ... #doctest: +NORMALIZE_WHITESPACE
    40.948 40.716

    """
    # Make mu and sigma unitless (user provides percentage)
    mu = mu / 100.0
    sigma = sigma / 100.0

    combs = {}
    for bolo in range(len(bolo_data)):
        sq = squid_ids[bolo]
        if sq not in combs:
            combs[sq] = []
        combs[sq].append((bolo_ids[bolo], bolo))

    tsout = 0.0 + bolo_data
    # How much to leak from one bolometer to its neighboring channels
    state = np.random.RandomState(seed)
    cross_amp = state.normal(mu, sigma, len(bolo_data))

    # Take one squid
    for sq1 in combs:
        # Take all bolometers in that SQUID
        for ch, i in combs[sq1]:
            # Take a second SQUID
            for sq2 in combs:
                if sq2 == sq1:
                    continue
                # Take all bolometers in that SQUID
                for ch2, i2 in combs[sq2]:
                    tsout[i] += cross_amp[i2] / squid_attenuation * tsout[i2]

    if new_array is not None:
        new_array[:] = tsout
    else:
        bolo_data[:] = tsout


def modify_beam_offsets(
    bolometer_xpos,
    bolometer_ypos,
    mu_diffpointing=10.0,
    sigma_diffpointing=5.0,
    seed=5847,
):
    """
    Modify the beam offsets (inject differential pointing between
    two pixel-pair bolometers). The model is the following:
    * Draw from a normal distribution G(mu, sigma)
      the magnitudes of the differential pointing rho.
    * Draw from a uniform distribution U(0, 2pi) the directions
      of the differential pointing theta.
    * Move the position of bottom bolometers as
      - x_top/bottom = \pm rho / 2 * cos(theta)
      - y_top/bottom = \pm rho / 2 * sin(theta)

    Parameters
    ----------
    bolometer_xpos : 1d array
        Initial beam centroid positions along the x axis. In radian.
    bolometer_ypos : 1d array
        Initial beam centroid positions along the y axis. In radian.
    mu_diffpointing : float, optional
        Mean of the distribution for differential pointing in arcsecond.
    sigma_diffpointing : float, optional
        Width of the distribution for differential pointing in arcsecond.
    seed : int, optional
        Random seed to control numbers generation.

    Returns
    ----------
    bolometer_xpos : 1d array
        Modified beam centroid positions along the x axis. In radian.
    bolometer_ypos : 1d array
        Modified beam centroid positions along the y axis. In radian.

    Examples
    ----------
    Inject differential pointing with magnitude 10 arcmin +- 5 arcmin.
    (completely irrealistic! just for test purposes).
    >>> bolometer_xpos = np.array([1., 1., -1., -1.]) * np.pi / 180.
    >>> bolometer_ypos = np.array([1., -1., -1., 1.]) * np.pi / 180.
    >>> x, y = modify_beam_offsets(bolometer_xpos, bolometer_ypos,
    ...     mu_diffpointing=600., sigma_diffpointing=300., seed=5847)
    >>> print(x * 180 / np.pi, y * 180 / np.pi) #doctest: +NORMALIZE_WHITESPACE
    [ 1.03802729  0.96197271 -1.02689298 -0.97310702]
    [ 0.92369044 -0.92369044 -1.10310547  1.10310547]

    """
    assert len(bolometer_xpos) == len(bolometer_ypos), ValueError(
        "x and y bolometer coordinates should have the same length"
    )
    state = np.random.RandomState(seed)

    npair = int(len(bolometer_xpos) / 2)

    # rho is defined by the user [xpos, ypos in radians - convert it!]
    rho = state.normal(
        mu_diffpointing * arcsecond2rad, sigma_diffpointing * arcsecond2rad, npair
    )

    # Angle are uniformly drawn from 0 and 360 degree
    theta = state.uniform(0, 360, npair) * deg2rad

    rhox = rho / 2.0 * np.cos(theta)
    rhoy = rho / 2.0 * np.sin(theta)

    bolometer_xpos[::2] += rhox
    bolometer_xpos[1::2] += -rhox
    bolometer_ypos[::2] += rhoy
    bolometer_ypos[1::2] += -rhoy

    return bolometer_xpos, bolometer_ypos


def inject_beam_ellipticity(
    sigma_gaussian,
    mu_beamellipticity,
    sigma_beamellipticity,
    nbolo,
    do_diffbeamellipticity=True,
    seed=54875,
):
    """
    Inject differential beam ellipticity (and beam ellipticity)
    Starting from the definition of the beam ellipticity:
        * eps = (sig1^2 - sig2^2)/(sig1^2 + sig2^2)
        * Assume ellipticity follows normal distribution
          centered on X% with Y% width
        * Assume sig1 = sig2 + d, where d is a small deviation
        * relate d and eps
    Furthermore, we assume the angle of the ellipses are taken
    from an uniform distribution btw -90 and 90 degrees.

    Parameters
    ----------
    sigma_gaussian : float
        Initial beam size (sigma). [radian]
    mu_beamellipticity : float
        Mean of the beam ellipticity distribution. [in percent]
        Between -100 and 100 percent.
    sigma_beamellipticity : float
        Width of the beam ellipticity distribution. [in percent]
    nbolo : int
        Number of bolometers for which beams have to be perturbed.
    do_diffbeamellipticity : bool, optional
        If False, bolometers in the same pair will have the same ellipticity.
        If True, add differential beam ellipticity between two bolometers.
        Default is True.
    seed : int, optional
        Seed for random number generation.

    Returns
    ----------
    sig_1 : 1d array of float
        Vector containing new values for semi-major axes for all bolometers
        in an elliptical Gaussian representation of the beam.
    sig_2 : 1d array of float
        Vector containing new values for semi-minor axes for all bolometers
        in an elliptical Gaussian representation of the beam.
    ellip_ang : 1d array of float
        Angles of rotation of beams around 0 for all bolometers.

    Examples
    ----------
    >>> from s4cmb.instrument import BeamModel
    >>> from s4cmb.instrument import FocalPlane
    >>> fp = FocalPlane(npair_per_squid=2, verbose=False)
    >>> bm = BeamModel(fp, verbose=False)

    Unperturbed beams (2 pairs of bolometers)
    >>> print(bm.sig_1)
    [ 0.00043235  0.00043235  0.00043235  0.00043235]

    Different pairs have different ellipticity
    (but not differential ellipticity within pairs)
    >>> sig_1, sig_2, ellip_ang = inject_beam_ellipticity(
    ...     bm.sig_1[0], 10, 5, 4, do_diffbeamellipticity=False)
    >>> print(sig_1)
    [ 0.0004715   0.0004715   0.00045276  0.00045276]

    Different pairs have different ellipticity and differential ellipticity
    within pairs
    >>> sig_1, sig_2, ellip_ang = inject_beam_ellipticity(
    ...     bm.sig_1[0], 10, 5, 4, do_diffbeamellipticity=True)
    >>> print(sig_1)
    [ 0.0004715   0.00045571  0.00045276  0.00046666]

    """
    state = np.random.RandomState(seed)

    eps = state.normal(
        mu_beamellipticity / 100.0, sigma_beamellipticity / 100.0, nbolo
    )

    if do_diffbeamellipticity:
        # d_plus = 2 * sigma_gaussian / eps * (1.0 + np.sqrt((1 - eps ** 2)))
        d_minus = 2 * sigma_gaussian / eps * (1.0 - np.sqrt((1 - eps ** 2)))
        d = d_minus
    else:
        # Bolometers in the same pair will have the same ellipticity.
        # This assumes that bolometers i and i+1 belong
        # to the same pair (default behaviour)
        eps = np.repeat(eps[::2], 2)
        # d_plus = 2 * sigma_gaussian / eps * (1.0 + np.sqrt((1 - eps ** 2)))
        d_minus = 2 * sigma_gaussian / eps * (1.0 - np.sqrt((1 - eps ** 2)))
        d = d_minus

    sig_1 = np.ones(nbolo) * sigma_gaussian + d / 2.0
    sig_2 = np.ones(nbolo) * sigma_gaussian - d / 2.0

    # Angle between ellipses
    ellip_ang = state.uniform(-90, 90, nbolo)

    return sig_1, sig_2, ellip_ang


def modify_pointing_parameters(values, errors):
    """
    This routine is not realistic at all for the moment!

    Parameters
    ----------
    values : 1d array
        Array containing values of the pointing parameters in degree.
    errors : 1d array
        Array containing values of the pointing parameter errors in degree.

    Returns
    ----------
    values_mod : 1d array
        Array containing modified values of the pointing parameters in degree.

    """
    values_mod = [p + err for p, err in zip(values, errors)]
    return values_mod


def step_function(nbolos, nsamples, mean=1, std=0.05, nbreaks=1, sign="same", seed=0):
    """
    Generate step functions for each bolometer from 1 to a values
    drawn from N(mean, std). The full timestream is broken into nbreaks
    corresponding to retuning of bolometers (gains are reset to one
    after a break). Return the gains for all bolometers (2 bolo in a
    pair will have the same break to avoid differential gain).

    Parameters
    ----------
    nbolos : int
        The number of bolometers
    nsamples : int
        Length of timestreams
    mean : float
        Mean of the distribution to sample drifts from. Default is 1.
    std : float
        Std of the distribution to sample drifts from.
        Default is 0.05 (that is \pm 5% with respect to a mean=1).
    nbreaks : int
        Number of break (number of retuning).
    sign : string
        If same, both detectors in the pair will have the same gain.
        If opposite, detectors will have gains with opposite sign
        (differential gain).

    Returns
    ----------
    gains : 2D array of size (nbolos, nsamples)
        The gains at each timestep for all bolometers

    Examples
    ----------
    >>> nbolos = 4
    >>> nsamples = 4
    >>> gains = step_function(nbolos, nsamples, nbreaks=1)
    >>> print(gains[0])
    [ 1.          1.          1.08820262  1.08820262]

    >>> gains = step_function(nbolos, nsamples, nbreaks=1, sign='opposite')
    >>> assert gains[0][0] == 2 - gains[1][0]
    """
    # Fix the seed
    state = np.random.RandomState(seed)

    # Initialise gains to ones
    gains = np.ones((nbolos, nsamples))

    # Length of each break
    length = nsamples // nbreaks
    sublength = nsamples // (2 * nbreaks)

    for pos in range(nbreaks):
        # 2 bolo in a pair will have the same break to avoid differential gain
        end_points = state.normal(mean, std, size=int(nbolos / 2))
        end_points = np.tile(end_points, (2, 1)).T.reshape((2 * len(end_points), 1))

        # Assign values
        shift = pos * length
        if pos * length < nsamples:
            gains[:, shift + sublength: shift + 2 * sublength] = end_points.reshape(
                (len(end_points), 1)
            )
        else:
            continue

        if sign == "opposite":
            gains[1] = 2 - gains[1]

    return gains


def step_function_gen(nsamples, mean=1, std=0.05, nbreaks=1, sign="same", seed=0):
    """
    Generator of step functions for each bolometer from 1 to a values
    drawn from N(mean, std). The full timestream is broken into nbreaks
    corresponding to retuning of bolometers (gains are reset to one
    after a break). Return the gains for 2 bolometers (2 bolo in a
    pair will have the same break to avoid differential gain).

    Particularly useful to save memory when simulating a lot of detectors.

    Parameters
    ----------
    nsamples : int
        Length of timestreams
    mean : float
        Mean of the distribution to sample drifts from. Default is 1.
    std : float
        Std of the distribution to sample drifts from.
        Default is 0.05 (that is \pm 5% with respect to a mean=1).
    nbreaks : int
        Number of break (number of retuning).
    sign : string
        If same, both detectors in the pair will have the same gain.
        If opposite, detectors will have gains with opposite sign
        (differential gain).

    Returns
    ----------
    gains : 2D array of size (2, nsamples)
        The gains at each timestep for 2 bolometers in a pair

    Examples
    ----------
    >>> nsamples = 4
    >>> gains_gen = step_function_gen(nsamples, nbreaks=1)
    >>> print(next(gains_gen)[0])
    [ 1.          1.          1.08820262  1.08820262]

    >>> gains_gen = step_function_gen(nsamples, nbreaks=1, sign='opposite')
    >>> g = next(gains_gen)
    >>> assert g[0][0] == 2 - g[1][0]
    """
    # Fix the seed
    state = np.random.RandomState(seed)

    # Initialise gains to ones
    gains = np.ones((2, nsamples))

    # Length of each break
    length = nsamples // nbreaks
    sublength = nsamples // (2 * nbreaks)

    while 1:
        for pos in range(nbreaks):
            # 2 bolo in a pair will have the same break
            # to avoid differential gain
            end_points = state.normal(mean, std, size=1)
            end_points = np.tile(end_points, (2, 1)).T.reshape((2, 1))

            # Assign values
            shift = pos * length
            if pos * length < nsamples:
                gains[
                    :, shift + sublength: shift + 2 * sublength
                ] = end_points.reshape((len(end_points), 1))
            else:
                continue

        if sign == "opposite":
            gains[1] = 2 - gains[1]

        yield gains


def linear_function(
    nbolos, nsamples, mean=1, std=0.05, nbreaks=1, sign="same", seed=0
):
    """
    Generate linear functions for each bolometer from 1 to a values
    drawn from N(mean, std). The full timestream is broken into nbreaks
    corresponding to retuning of bolometers (gains are reset to one
    after a break). Return the gains for all bolometers (2 bolo in a
    pair will have the same break to avoid differential gain)

    Parameters
    ----------
    nbolos : int
        The number of bolometers
    nsamples : int
        Length of timestreams
    mean : float
        Mean of the distribution to sample drifts from. Default is 1.
    std : float
        Std of the distribution to sample drifts from.
        Default is 0.05 (that is \pm 5% with respect to a mean=1).
    nbreaks : int
        Number of break (number of retuning).
    sign : string
        If same, both detectors in the pair will have the same gain.
        If opposite, detectors will have gains with opposite sign
        (differential gain).

    Returns
    ----------
    gains : 2D array of size (nbolos, nsamples)
        The gains at each timestep for all bolometers

    Examples
    ----------
    >>> nbolos = 4
    >>> nsamples = 4
    >>> gains = linear_function(nbolos, nsamples, nbreaks=1)
    >>> print(gains[0])
    [ 1.          1.02940087  1.05880174  1.08820262]

    >>> gains = linear_function(nbolos, nsamples, nbreaks=1, sign='opposite')
    >>> assert gains[0][0] == 2 - gains[1][0]
    """
    # Fix the seed
    state = np.random.RandomState(seed)

    # Initialise gains to ones
    gains = np.ones((nbolos, nsamples))

    # Length of each break
    length = nsamples // nbreaks

    for pos in range(nbreaks):
        # 2 bolo in a pair will have the same break to avoid differential gain
        end_points = state.normal(mean, std, size=int(nbolos / 2))
        end_points = np.tile(end_points, (2, 1)).T.reshape((2 * len(end_points), 1))

        # Assign values
        shift = pos * length
        if pos * length < nsamples:
            gains[:, shift: shift + length] = np.array(
                [
                    np.interp(
                        range(shift, shift + length),
                        [shift, shift + length - 1],
                        [1, end[0]],
                    )
                    for end in end_points
                ]
            )
        else:
            continue

    if sign == "opposite":
        gains[1] = 2 - gains[1]

    return gains


def linear_function_gen(nsamples, mean=1, std=0.05, nbreaks=1, sign="same", seed=0):
    """
    Generator returning linear functions for each bolometer from 1 to a values
    drawn from N(mean, std). The full timestream is broken into nbreaks
    corresponding to retuning of bolometers (gains are reset to one
    after a break). Yield the gains for 2 bolometers (2 bolo in a
    pair will have the same break to avoid differential gain).

    Particularly useful to save memory when simulating a lot of detectors.

    Parameters
    ----------
    nsamples : int
        Length of timestreams
    mean : float
        Mean of the distribution to sample drifts from. Default is 1.
    std : float
        Std of the distribution to sample drifts from.
        Default is 0.05 (that is \pm 5% with respect to a mean=1).
    nbreaks : int
        Number of break (number of retuning).
    sign : string
        If same, both detectors in the pair will have the same gain.
        If opposite, detectors will have gains with opposite sign
        (differential gain).

    Returns
    ----------
    gains : 2D array of size (2, nsamples)
        The gains at each timestep for 2 bolometers in a pair

    Examples
    ----------
    >>> nsamples = 4
    >>> gains_gen = linear_function_gen(nsamples, nbreaks=1)
    >>> print(next(gains_gen)[0])
    [ 1.          1.02940087  1.05880174  1.08820262]

    >>> gains_gen = linear_function_gen(nsamples, nbreaks=1, sign='opposite')
    >>> g = next(gains_gen)
    >>> assert g[0][0] == 2 - g[1][0]

    """
    # Fix the seed
    state = np.random.RandomState(seed)

    # Initialise gains to ones
    gains = np.ones((2, nsamples))

    # Length of each break
    length = nsamples // nbreaks

    while 1:
        for pos in range(nbreaks):
            # 2 bolo in a pair will have the same
            # break to avoid differential gain
            end_points = state.normal(mean, std, size=1)
            end_points = np.tile(end_points, (2, 1)).T.reshape((2, 1))

            # Assign values
            shift = pos * length
            if pos * length < nsamples:
                gains[:, shift: shift + length] = np.array(
                    [
                        np.interp(
                            range(shift, shift + length),
                            [shift, shift + length - 1],
                            [1, end[0]],
                        )
                        for end in end_points
                    ]
                )
            else:
                continue

        if sign == "opposite":
            gains[1] = 2 - gains[1]

        yield gains


def get_kernel_coefficients(
    beamprm, pairlist, kernel_type="diff", nx=128, pix_size=None, basis="circular"
):
    """
    Generate beam expansion coefficients.

    Here we compute the coefficients of the operator K to go from the observed
    temperature map to the T->P and P->T as well as T->T, P->P leakage. It is
    based on the sum and difference beam maps derivative expansion onto
    a basis beam:

    diffbeam = K- conv basis_beam,
    sumbeam  = K+ conv basis_beam,

    with

    K+/- = a.D,

    where D is a vector containing beam map and its derivatives:

    D = (I, dI/dt, dI/dp, d2I/dpdt, d2I/dt2, d2I/d2p),

    and a are coefficients. This assumes the polarized beams are the same for
    intensity and linear polarization.

    Parameters
    ----------
    beamprm : beam_model instance
        Instance of beam_model.
    pairlist : list of list
        List containing the indices of bolometers grouped by pair
        [[0, 1], [2, 3], ...]
    kernel_type : string
        Specify if the kernel needed ar for a sum or diff beams or both.
        Supported options: 'sum','diff' or both.
    nx : int, optional
        Number of pixels per row/column (in pixel) to construct the beam maps.
        You want nx * pix_size to be a large number compared to the beam size
        to make sure to incorporate all beam features.
    pix_size : float, optional
        Pixel size in radian. If None, set automatically to be 1/7th of the
        beam size.
    basis : string
        Specify the basis onto whish developing the sum or diff beam. The
        default assumes a Gaussian circular beam with FWHM specified in the beam
        model instance. The sum beam can be used instead.
        Supported options: 'circular' (default),'sumbeam'

    Returns
    ----------
    out : array of array
        Array containing beam kernel coefficients
        (one for each T derivative) for all pairs.

    Examples
    ----------
    >>> from s4cmb.instrument import FocalPlane
    >>> from s4cmb.instrument import BeamModel
    >>> fp = FocalPlane(verbose=False)
    >>> pairlist = np.reshape(fp.bolo_index_in_fp, (fp.npair, 2))

    If the beams are identical, the coefficients will be zeros
    >>> bm = BeamModel(fp, verbose=False)
    >>> coeffs = get_kernel_coefficients(bm, pairlist, nx=32)
    >>> np.allclose(coeffs[0], np.zeros_like(coeffs[0]))
    True

    Let's perturb the beam now by adding beam ellipticity.
    Beam ellipticity involves a Stokes parameter and its second derivatives,
    but not the first derivatives (2nd and 3rd coefficients are zeros).
    >>> bm = BeamModel(fp, verbose=False)
    >>> bm.sig_1, bm.sig_2, bm.ellip_ang = inject_beam_ellipticity(
    ...     bm.sig_1[0], 10, 5, fp.nbolometer, do_diffbeamellipticity=True)
    >>> coeffs = get_kernel_coefficients(bm, pairlist, nx=32)

    Let's perturb the beam now by adding differential pointing.
    Differential pointing mainly the first derivatives of a Stokes parameter,
    but not the temperature itself and its second derivatives
    (0th, 4th, 5th, and 6th coefficients are zeros).
    >>> bm = BeamModel(fp, verbose=False)
    >>> bm.xpos, bm.ypos = modify_beam_offsets(bm.xpos, bm.ypos,
    ...     mu_diffpointing=600., sigma_diffpointing=300., seed=5847)
    >>> coeffs = get_kernel_coefficients(bm, pairlist, nx=32)

    """
    assert basis in ["circular", "sumbeam"]
    assert kernel_type in ["sum", "diff", "sumdiff", "diffsum"]
    if pix_size is None:
        # Go from sigma to FWHM
        mean_FWHM_x = np.mean(beamprm.sig_1) / np.sqrt(8 * np.log(2))
        mean_FWHM_y = np.mean(beamprm.sig_2) / np.sqrt(8 * np.log(2))

        # 1/7th of the beam size.
        pix_size = (mean_FWHM_x + mean_FWHM_y) / 2.0 / 7.0

    if basis == "circular":
        # Creates the a circular beam used to develop the sum/diff beam
        fwhm = beamprm.fwhm * np.pi / 180 / 60
        xy2f = coordinates_on_grid(pix_size=pix_size, nx=nx)
        basis_beam = gauss2d(
            xy2f,
            0,
            0,
            1,
            fwhm / np.sqrt(8 * np.log(2)),
            fwhm / np.sqrt(8 * np.log(2)),
            0,
        ).reshape((nx, nx))

    out_diff = []
    out_sum = []

    for ct, cb in pairlist:
        summap, diffmap = construct_beammap(beamprm, ct, cb, nx, pix_size)
        if summap is not None and diffmap is not None:
            if basis != "circular":
                basis_beam = summap
            if "diff" in kernel_type:
                kernel_diff = split_deriv(basis_beam, diffmap, pix_size)
            if "sum" in kernel_type:
                kernel_sum = split_deriv(basis_beam, summap, pix_size)
        else:
            kernel_sum = None
        if "diff" in kernel_type:
            out_diff.append(kernel_diff)
        if "sum" in kernel_type:
            out_sum.append(kernel_sum)
    if ("diff" in kernel_type) and ("sum" in kernel_type):
        return np.array(out_diff), np.array(out_sum)
    else:
        if "diff" == kernel_type:
            return np.array(out_diff)
        if "sum" == kernel_type:
            return np.array(out_sum)


def split_deriv(sumbeam, diffbeam, pix_size):
    """
    Try to split a convolution kernel diffbeam (B-) in to two convolutions,
    one low order derivatives (K), and the other a base kernel (B+):

    B- = K conv B+

    where K = a.D, with D is a vector containing Stokes parameter map
    and its derivatives. In the case of e.g. temperature:

    D = (I, dI/dt, dI/dp, d2I/dpdt, d2I/dt2, d2I/d2p),

    and a are coefficients.

    Parameters
    ----------
    sumbeam : 2d array
        Sum beam map
    diffbeam : 2d array
        Difference beam map
    pix_size : float
        Size of a pixel [radian]

    Returns
    ----------
    x : 1d array of 6 elements
        Vector containing the beam kernel coefficients a
        (one for each derivative of the beam map).

    Examples
    ----------
    >>> from s4cmb.instrument import FocalPlane
    >>> from s4cmb.instrument import BeamModel
    >>> fp = FocalPlane(verbose=False)
    >>> bm = BeamModel(fp, verbose=False)
    >>> pix_size = 0.5 / 60. * np.pi / 180.
    >>> summap, diffmap = construct_beammap(bm, 0, 1, nx=32, pix_size=pix_size)
    >>> K = split_deriv(summap, diffmap, pix_size)

    """
    ds = derivs(sumbeam, pix_size)

    A = ds.reshape((6, sumbeam.shape[0] * sumbeam.shape[1])).T
    b = diffbeam.reshape(diffbeam.shape[0] * diffbeam.shape[1])

    x, residues, rank, s = np.linalg.lstsq(A, b, rcond=-1)

    return x


def xderiv(m, pix_size):
    """
    Perform the derivative of a 2d map with respect to x coordinate.

    Parameters
    ----------
    m : 2d array
        Input 2d map
    pix_size : float
        Size of a pixel [radian]

    Returns
    ----------
    -v : 2d array
        Derivative of the input beam map with respect to x coordinate.

    Examples
    ----------
    >>> m = np.ones((10, 10)) * np.arange(10)
    >>> dmx = xderiv(m, pix_size=0.5/60.*np.pi/180.)
    """
    kernel = np.array([[1, 0, -1]]) / (2.0 * pix_size)
    v = signal.convolve2d(m, kernel, mode="same")
    return -v


def yderiv(m, pix_size):
    """
    Perform the derivative of a 2d map with respect to y coordinate.

    Parameters
    ----------
    m : 2d array
        Input 2d map
    pix_size : float
        Size of a pixel [radian]

    Returns
    ----------
    -v : 2d array
        Derivative of the input beam map with respect to y coordinate.

    Examples
    ----------
    >>> m = np.ones((10, 10)) * np.arange(10)
    >>> dmy = yderiv(m.T, pix_size=0.5/60.*np.pi/180.)
    """
    kernel = np.array([[1, 0, -1]]).T / (2.0 * pix_size)
    v = signal.convolve2d(m, kernel, mode="same")
    return -v


def derivs(m, pix_size):
    """
    Compute full 1st and 2nd derivatives of a 2d map (flat sky).

    Parameters
    ----------
    m : 2d array
        Input 2d map
    pix_size : float
        Size of a pixel [radian]

    Returns
    ----------
    v : array of 2d array
        Derivatives of the input beam map with respect to x and y coordinates.

    Examples
    ----------
    >>> m = np.ones((10, 10)) * np.arange(10)
    >>> m = m.T * np.arange(10)
    >>> m00, m10, m01, m11, m20, m02 = derivs(m, pix_size=0.5/60.*np.pi/180.)
    """
    m00 = m
    m10 = xderiv(m, pix_size)
    m01 = yderiv(m, pix_size)
    m11 = yderiv(m10, pix_size)
    m20 = xderiv(m10, pix_size)
    m02 = yderiv(m01, pix_size)

    return np.array((m00, m10, m01, m11, m20, m02))


def fixspin(K, spin):
    """
    Selects derivatives of only a given spin (by default we do not compute
    all spins if this is unecessary) and rearrange for later use.

    Parameters
    ----------
    K : 1d array
        Kernel coefficients of the beam map decomposition
    spins : string
        Which type of derivatives you want to use for you leakage.
        Example: beam ellipticity would use 0 and 2, while differential
        pointing would use just 1.

    Returns
    ----------
    es : 1d array
        Selected kernel coefficients of the beam map decomposition

    Examples
    ----------
    >>> from s4cmb.instrument import FocalPlane
    >>> from s4cmb.instrument import BeamModel
    >>> fp = FocalPlane(verbose=False)
    >>> bm = BeamModel(fp, verbose=False)
    >>> pix_size = 0.5 / 60. * np.pi / 180.
    >>> summap, diffmap = construct_beammap(bm, 0, 1, nx=32, pix_size=pix_size)
    >>> K = split_deriv(summap, diffmap, pix_size)
    >>> SK = fixspin(K, spin='012')

    """
    d00, d10, d01, d11, d20, d02 = K
    x = (d20 + d02) * 0.5
    y = (d20 - d02) * 0.5

    e00 = 0.0
    e10 = 0.0
    e01 = 0.0
    e11 = 0.0
    e20 = 0.0
    e02 = 0.0

    if "0" in spin:
        e00 += d00
        e10 += 0.0
        e01 += 0.0
        e11 += 0.0
        e20 += x
        e02 += x
    if "1" in spin:
        e00 += 0.0
        e10 += d10
        e01 += d01
        e11 += 0.0
        e20 += 0.0
        e02 += 0.0
    if "2" in spin:
        e00 += 0.0
        e10 += 0.0
        e01 += 0.0
        e11 += d11
        e20 += y
        e02 += -y

    return np.array((e00, e10, e01, e11, e20, e02))


def rotate_deriv(K, theta):
    """
    Rotate a vector of derivative kernel coefficients K by an angle theta.

    Parameters
    ----------
    K : 1d array
        Kernel coefficients of the beam map decomposition
    theta : float
        Orientation of the pixel [radian]

    Returns
    ----------
    es : 1d array
        Rotated kernel coefficients of the beam map decomposition.

    Examples
    ----------
    >>> from s4cmb.instrument import FocalPlane
    >>> from s4cmb.instrument import BeamModel
    >>> fp = FocalPlane(verbose=False)
    >>> bm = BeamModel(fp, verbose=False)
    >>> pix_size = 0.5 / 60. * np.pi / 180.
    >>> summap, diffmap = construct_beammap(bm, 0, 1, nx=32, pix_size=pix_size)
    >>> K = split_deriv(summap, diffmap, pix_size)
    >>> RK = rotate_deriv(K, theta=np.pi)

    """

    d00, d10, d01, d11, d20, d02 = K

    c = np.cos(theta)
    s = np.sin(theta)

    e00 = np.ones_like(theta) * d00

    e10 = c * d10 - s * d01
    e01 = s * d10 + c * d01

    e20 = c * c * d20 - 2.0 * c * s * d11 + s * s * d02
    e02 = s * s * d20 + 2.0 * c * s * d11 + c * c * d02
    e11 = (c * c - s * s) * d11 + c * s * (d20 - d02)

    es = np.array((e00, e10, e01, e11, e20, e02))

    return es


def fix_kernel_type(K, kernel_type="diffnomonopole"):
    if "sum" in kernel_type:
        # subtracts monopole part already included
        # in the input timestream
        K[:, 0] -= 1
    if "nomonopole" in kernel_type:
        # Force monopole to 0 when using leakage beam
        K[:, 0] = 0.0
    # Flips the sign of dtheta terms.
    K[:, 1] *= -1
    K[:, 3] *= -1
    return K


def waferts_add_diffbeam(
    waferts,
    point_matrix,
    beam_orientation,
    intensity_derivatives,
    diffbeam_kernels,
    pairlist,
    spins="012",
    kernel_type="diff_nomonopole",
    pol_derivatives=None,
    pol_angle=None,
    spins_pol="012",
    Usign=1,
):
    """
    Modify timestreams by injecting T->P leakage from beam mismatch.
    Note that timestreams are modified on-the-fly.

    Parameters
    ----------
    waferts : ndarray
        Array containing timestreams
    point_matrix : ndarray
        Array containing pointing information
    beam_orientation : ndarray
        Array containing intrinsic orientation of the pixels
    intensity_derivatives : ndarray
        Containing map of the intensity and its derivatives (6 maps). They need
        to include only the pixels in the field defined by the scan (see obspix
        field of TOD class).
    diffbeam_kernels : 1d array of 6 elements
        Array containing beam kernel coefficients
        (one for each T derivative) for all pairs.
    pairlist : list of list
        List containing the indices of bolometers grouped by pair
        [[0, 1], [2, 3], ...].
    spins : string
        Which type of derivatives you want to use for you leakage.
        Example: beam ellipticity would use 0 and 2, while differential
        pointing would use just 1.
    kernel_type : string
        String specifying if sum or diff beam leakage components
        has to be injected. Add _nomonopole if you want not to include
        monopole leakage component and force it to be 0.
    pol_derivatives : list of ndarray
        Containing map of the Q and U parameters and their
        derivatives (2x6 maps: [[Q,dQ/dx...],[U,dU/dx...]])
    pol_angle: ndarray
        Array containing the value of the polarization angle
    spins_pol: string
        as the spins argument. Allows to use different leakage
        derivatives for the polarization part of the beam leakage. By
        default spins=spins_pol
    Usign: int
        defines sign convention to match the IAU convention in case
        the flat map projection is used.

    Examples
    ----------
    >>> from s4cmb.tod import load_fake_instrument
    >>> from s4cmb.tod import TimeOrderedDataPairDiff
    >>> inst, scan, sky_in = load_fake_instrument(compute_derivatives=True)

    Make beam ellipticity
    >>> sig_1, sig_2, ellip_ang = inject_beam_ellipticity(
    ...     inst.beam_model.sig_1[0], 50, 10,
    ...     inst.focal_plane.nbolometer, do_diffbeamellipticity=True)
    >>> inst.beam_model.sig_1 = sig_1
    >>> inst.beam_model.sig_2 = sig_2
    >>> inst.beam_model.ellip_ang = ellip_ang

    Get TOD
    >>> tod = TimeOrderedDataPairDiff(inst, scan, sky_in, CESnumber=0)
    >>> d = np.array([tod.map2tod(det) for det in range(2 * tod.npair)])

    Get beam kernels
    >>> pairlist = np.reshape(
    ...     inst.focal_plane.bolo_index_in_fp, (inst.focal_plane.npair, 2))
    >>> K = get_kernel_coefficients(
    ...     inst.beam_model, pairlist, nx=128, pix_size=None)

    Get intensity derivatives and orientation of pixels
    >>> intensity_derivatives = np.array(
    ...     [sky_in.I, sky_in.dIdt, sky_in.dIdp,
    ...      sky_in.d2Idpdt, sky_in.d2Id2t, sky_in.d2Id2p])
    >>> beam_orientation = np.array(
    ...     [tod.pol_angs[ch] - (
    ...         90.0 - tod.intrinsic_polangle[2*ch]) * np.pi / 180. -
    ...      2*tod.hwpangle for ch in range(inst.focal_plane.npair)])

    Inject spurious signal
    >>> waferts_add_diffbeam(
    ...     d, tod.point_matrix, beam_orientation,
    ...     intensity_derivatives, K,
    ...     pairlist, spins='012')

    """
    assert kernel_type in ["sum", "diff", "sum_nomonopole", "diff_nomonopole"]

    # Just to preserve it
    K = diffbeam_kernels + 0.0

    K = fix_kernel_type(K, kernel_type)

    if "sum" in kernel_type:
        # adds I->I in sum timestream and P->P in diff timestream
        bottom_sign = -1
    else:
        bottom_sign = 1

    for i in range(len(K)):
        K[i] = fixspin(K[i], spins)

    diffbeamleak = np.zeros((int(waferts.shape[0] / 2), waferts.shape[1]))

    # adds I leakage. if kernel_type is "diff" adds I->P leakage, otherwise
    # adds I->I leakage
    diffbeam_map2tod(
        diffbeamleak, intensity_derivatives, point_matrix, beam_orientation, K
    )

    waferts[::2] += diffbeamleak
    waferts[1::2] -= bottom_sign * diffbeamleak

    # adds P leakage. If kernel_type is "diff" adds P->I leakage, otherwise
    # adds P->P leakage
    if (pol_derivatives is not None) and (pol_angle is not None):
        if (spins_pol != spins) and (pol_derivatives is not None):
            # initialize a new copy of polarization kernels withot modifying
            # the input
            K_pol = fix_kernel_type(diffbeam_kernels + 0.0, kernel_type)
            K_pol[i] = fixspin(K_pol[i], spins_pol)
        else:
            K_pol = K
        diffbeamleak[:, :] = 0.0
        diffbeam_map2tod(
            diffbeamleak,
            pol_derivatives,
            point_matrix,
            beam_orientation,
            K_pol,
            pol_angle,
            Usign,
        )

        waferts[::2] += diffbeamleak
        waferts[1::2] += bottom_sign * diffbeamleak


def diffbeam_map2tod(
    out,
    signal_derivatives,
    point_matrix,
    beam_orientation,
    diffbeam_kernels,
    pol_angle=None,
    Usign=1,
):
    """
    Scan the leakage maps to generate polarisation timestream for channel ch,
    using the differential beam model kernel.


    Note for future:
     * check flat sky!

    Parameters
    ----------
    out : ndarray
        Will contain the spurious signal (npair, nsamples)
    signal_derivatives : ndarray
        Containing map of the intensity and its derivatives (6 maps)
        or a list of Q and U derivatives (2x6 maps)
    point_matrix : ndarray
        Array containing pointing information
    beam_orientation : ndarray
        Array containing intrinsic orientation of the pixels
    diffbeam_kernels : 1d array of 6 elements
        Array containing beam kernel coefficients
        (one for each T derivative) for all pairs.
    pol_angle: ndarray
        Array containing polarization angle information.
    Usign: int
        defines sign convention to match the IAU convention in case
        the flat map projection is used.
    """
    npix, nt = out.shape
    assert point_matrix.shape == out.shape
    assert diffbeam_kernels.shape[0] == npix
    if pol_angle is None:
        assert len(signal_derivatives) == diffbeam_kernels.shape[1]
    else:
        assert len(signal_derivatives) == 2
        assert len(signal_derivatives[0]) == diffbeam_kernels.shape[1]
    assert beam_orientation.shape == (npix, nt)

    for ipix in range(npix):
        i1d = point_matrix[ipix, :]
        okpointing = i1d != -1
        io = i1d[okpointing]
        k = rotate_deriv(diffbeam_kernels[ipix], -beam_orientation[ipix, okpointing])
        if pol_angle is not None:
            # compute cosine and sine polarization angles to save time
            psi = pol_angle[ipix, :]
            c2pol_ang = np.cos(2 * psi[okpointing])
            s2pol_ang = np.sin(2 * psi[okpointing])
            if signal_derivatives is not None:
                for coeff in range(len(k)):
                    # Q derivatives
                    der0 = signal_derivatives[0][coeff, io]
                    out[ipix, okpointing] += der0 * k[coeff] * c2pol_ang
                    # U derivatives
                    der1 = signal_derivatives[1][coeff, io]
                    out[ipix, okpointing] += der1 * k[coeff] * s2pol_ang * Usign
        # Add the leakage on-the-fly
        else:
            for coeff in range(len(k)):
                # T derivatives
                out[ipix, okpointing] += signal_derivatives[coeff, io] * k[coeff]


if __name__ == "__main__":
    import doctest

    if np.__version__ >= "1.14.0":
        np.set_printoptions(legacy="1.13")
    doctest.testmod()
