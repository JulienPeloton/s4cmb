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
Module including diverse tools to manipulate alms and maps .
"""
import numpy as np
import healpy as hp


def get_healpix_ring_pixel_layout(nside, th_idx):
    """Healpix ring layout.

    From 'get_pixel_layout' subroutine in healpix f90 package.

    Author: Julien Carron (j.carron@sussex.ac.uk)

    Parameters
    ----------
    nside : int
        Healpix nside parameter.
    th_idx : int
        Ring index (0 <= th_idx < 4*nside - 1).

    Returns
    -------
    startpix : int
        Starting pixel number.
    nphi : int
        Number of pixels in the ring.
    kphi0 : int
        Starting pixel number in the ring.
    cth : float
        Cosine of the polar angle.
    sth : float
        Sine of the polar angle.
    """
    ith = th_idx + 1
    nrings = 2 * nside
    assert 1 <= ith <= 4 * nside - 1, (ith, nrings)
    if ith > nrings:
        startpix, nphi, kphi0, cth, sth = get_healpix_ring_pixel_layout(
            nside, ith - 2 * (ith - nrings) - 1
        )
        return 12 * nside ** 2 - startpix - nphi, nphi, kphi0, -cth, sth
    dth1 = 1.0 / 3.0 / nside ** 2
    dth2 = 2.0 / 3.0 / nside
    dst1 = 1.0 / (np.sqrt(6.0) * nside)
    if ith < nside:  # polar cap (north)
        cth = 1.0 - ith ** 2 * dth1
        nphi = 4 * ith
        kphi0 = 1
        sth = np.sin(2.0 * np.arcsin(ith * dst1))
        startpix = 2 * ith * (ith - 1)
    else:
        cth = (2 * nside - ith) * dth2
        nphi = 4 * nside
        kphi0 = (ith + 1 - nside) % 2
        sth = np.sqrt((1.0 - cth) * (1.0 + cth))
        startpix = 2 * nside * (nside - 1) + (ith - nside) * int(nphi)
    return startpix, nphi, kphi0, cth, sth


def get_alpha_raise(s, lmax):
    """Response coefficient of spin-s spherical harmonic to spin raising
    operator.

    Author: Julien Carron (j.carron@sussex.ac.uk)

    Parameters
    ----------
    s : int
        Input spin of the spherical harmonic.
    lmax : int
        Maximum multipole moment.

    Returns
    -------
    ret : np.ndarray
        Response coefficient of spin-s spherical harmonic to spin raising
        operator.

    Notes
    -----
    The response coefficient is defined as:
        alpha(s, l) = sqrt((l - s) * (l + s + 1))
    where l is the multipole moment.
    The response coefficient is zero for l < |s|.
    """
    ret = np.zeros(lmax + 1, dtype=float)
    ret[abs(s):] = np.sqrt(
        np.arange(
            abs(s) - s,
            lmax - s + 1
        ) * np.arange(
            abs(s) + s + 1,
            lmax + s + 2
        )
    )
    return ret


def get_alpha_lower(s, lmax):
    """Response coefficient of spin-s spherical harmonic to spin lowering
    operator.

    Author: Julien Carron (j.carron@sussex.ac.uk)

    Parameters
    ----------
    s : int
        Input spin of the spherical harmonic.
    lmax : int
        Maximum multipole moment.

    Returns
    -------
    ret : np.ndarray
        Response coefficient of spin-s spherical harmonic to spin lowering
        operator.

    Notes
    -----
    The response coefficient is defined as:
        alpha(s, l) = sqrt((l + s) * (l - s + 1))
    where l is the multipole moment.
    The response coefficient is zero for l < |s|.
    """
    ret = np.zeros(lmax + 1, dtype=float)
    ret[abs(s):] = -np.sqrt(
        np.arange(
            s + abs(s),
            lmax + s + 1
        ) * np.arange(
            abs(s) - s + 1,
            lmax - s + 2
        )
    )
    return ret


def alm2map_spin_der1(gclm, nside, spin, zbounds=(-1.0, 1.0), ret_slice=None):
    """Returns spin-s transform '_{s}d' of alm,
    together with d/dtheta _{s}d and 1/sin tht d/dphi _{s}d.

    This crude version has three calls to spin-weight harmonics alm2map_spin.

    Author: Julien Carron (j.carron@sussex.ac.uk)

    Parameters
    ----------
    gclm : list[np.ndarray] or np.ndarray
        List of two arrays containing the "gradient" and "curl" parts of the
        spherical harmonic coefficients.
        gclm[0] is the part defined by -(alm_{|s|} + (-1)^s alm_{-|s|})/2 and
        gclm[1] is part given by -(alm_{|s|} + (-1)^s alm_{-|s|})/(2i),
        and the last dimension is the healpix ordering
        of the alms.
    nside : int
        Healpix nside parameter.
    spin : int
        Input spin of the fields given in gclm, must be >= 0.
    zbounds : tuple of ndarray
        Bounds of the polar angle (cosine) for the output map, must be of
        shape (2,), default is (-1.0, 1.0).
    ret_slice : slice, optional
        Slice to return a subset of the output map, the computation is however
        done on the entire healpix grid, default is None.

    Returns
    -------
    _sd : np.ndarray
        Spin-s maps of the input spherical harmonic, shape (2, 12 * nside ** 2)
        if ret_slide is None, otherwise (2, len(ret_slide)).
        Note that (_sd[0] + 1j*_sd[1]) is the spin (+s) field,
        and (_sd[0] - 1j*_sd[1]) is the spin (-s) field.
    d_dth : np.ndarray
        First partial derivative of the spin-s alms with respect to theta,
        shape (2, 12 * nside ** 2) if ret_slide is None,
        otherwise (2, len(ret_slide)).
    d_dphi_sin0 : np.ndarray
        First partial derivative of the spin-s alms with respect to phi with
        factor 1/(sin (theta)), shape (2, 12 * nside ** 2)
        if ret_slide is None, otherwise (2, len(ret_slide)).
    """

    assert spin >= 0, spin
    assert hp.Alm.getlmax(gclm[0].size) == hp.Alm.getlmax(gclm[1].size)
    lmax = hp.Alm.getlmax(gclm[0].size)
    zbounds = np.sort(np.array(zbounds))
    # shape (2, 12 * nside ** 2),

    if spin == 0:
        # If spin == 0, we can use the alm2map function on the input alms
        _sd, d_dth, d_dphi_sin0 = hp.alm2map_der1(
            -(gclm[0] + 1j*gclm[1]),
            nside,
            lmax)

        return (
            np.vstack([_sd.real, _sd.imag]),
            np.vstack([d_dth.real, d_dth.imag]),
            np.vstack([d_dphi_sin0.real, d_dphi_sin0.imag])
        )

    # first entry = real part, second entry imaginary part.

    # Computing the maps from the input alms,
    # where _sd[0] + 1j*_sd[1] is the spin (+s) fiedl,
    # and _sd[0] - 1j*_sd[1] is the spin (-s) field.
    _sd = -np.array([hp.alm2map(alms, nside) for alms in gclm])

    # First computing the spin-lowering operator on the input alms
    _gclm = [
            hp.almxfl(gclm[0], get_alpha_lower(spin, lmax)),
            hp.almxfl(gclm[1], get_alpha_lower(spin, lmax)),
        ]
    if spin > 1:
        _sm1d = np.array(hp.alm2map_spin(_gclm, nside, spin - 1, lmax))
    else:
        _sm1d = -np.array([hp.alm2map(alms, nside) for alms in gclm])

    # Then computing the spin-raising operator on the input alms
    _gclm = [
        hp.almxfl(gclm[0], get_alpha_raise(spin, lmax)),
        hp.almxfl(gclm[1], get_alpha_raise(spin, lmax)),
    ]
    _sp1d = np.array(hp.alm2map_spin(_gclm, nside, spin + 1, lmax))

    # Retrieving the application of the partial derivative over theta from
    # the spin-lowering and spin-raising operators
    d_dth = -0.5 * (_sp1d + _sm1d)

    # Retrieving the application of the partial derivative over phi
    # (with factor 1/sin(theta)) from the spin-lowering and spin-raising
    # operators
    d_dphi_sin0 = 0.5 * np.array([-_sp1d[1] + _sm1d[1], _sp1d[0] - _sm1d[0]])
    for iring in range(4 * nside - 1):
        startpix, nphi, kphi0, cth, sth = get_healpix_ring_pixel_layout(
            nside,
            iring
        )
        if zbounds[0] <= cth <= zbounds[1]:
            slic = slice(startpix, startpix + nphi)
            d_dphi_sin0[1, slic] -= spin * (cth / sth) * _sd[0, slic]
            d_dphi_sin0[0, slic] += spin * (cth / sth) * _sd[1, slic]
    if ret_slice is not None:
        return (
            _sd[:, ret_slice],
            d_dth[:, ret_slice],
            d_dphi_sin0[:, ret_slice])

    return _sd, d_dth, d_dphi_sin0


def alm2map_spin_der2(gclm, nside, spin, zbounds=(-1.0, 1.0), ret_slice=None):
    """Returns spin-s transform '_{s}d' of alm, together with
        * d/dtheta _{s}d
        * 1/sin(theta) d/dphi _{s}d
        * d^2/dtheta^2 _{s}d
        * 1/sin(theta) d/dphi d/dtheta _{s}d
        * 1/(sin(theta)^2) d^2/dphi^2 _{s}d
    given as functions of the polar angle (cosine) theta.

    Parameters
    ----------
    gclm : list[np.ndarray] or np.ndarray
        List of two arrays containing the "gradient" and "curl" parts of the
        spherical harmonic coefficients.
        gclm[0] is the part defined by -(alm_{|s|} + (-1)^s alm_{-|s|})/2 and
        gclm[1] is part given by -(alm_{|s|} + (-1)^s alm_{-|s|})/(2i),
        and the last dimension is the healpix ordering
        of the alms.
    nside : int
        Healpix nside parameter.
    spin : int
        Input spin of the fields given in gclm, must be >= 0.
    zbounds : tuple of ndarray
        Bounds of the polar angle (cosine) for the output map, must be of
        shape (2,), default is (-1.0, 1.0).
    ret_slice : slice, optional
        Slice to return a subset of the output map, the computation is however
        done on the entire healpix grid, default is None.

    Returns
    -------
    _sd : np.ndarray
        Spin-s maps of the input spherical harmonic, shape (2, 12 * nside ** 2)
        if ret_slide is None, otherwise (2, len(ret_slide)).
        Note that (_sd[0] + 1j*_sd[1]) is the spin (+s) field,
        and (_sd[0] - 1j*_sd[1]) is the spin (-s) field.
    d_dth : np.ndarray
        First partial derivative of the spin-s alms with respect to theta,
        shape (2, 12 * nside ** 2) if ret_slide is None,
        otherwise (2, len(ret_slice)).
    d_dphi_sin0 : np.ndarray
        First partial derivative of the spin-s alms with respect to phi with
        factor 1/(sin theta), shape (2, 12 * nside ** 2) if ret_slide is None,
        otherwise (2, len(ret_slice)).
    d2_dth2 : np.ndarray
        Second partial derivative of the spin-s alms with respect to theta,
        shape (2, 12 * nside ** 2) if ret_slide is None,
        otherwise (2, len(ret_slice)).
    d_phi_d_th : np.ndarray
        Succession of partial derivatives of the spin-s alms with respect
        to phi and theta with factor 1/(sin theta),
        shape (2, 12 * nside ** 2)
        if ret_slide is None, otherwise (2, len(ret_slice)).
    d2_dphi2_sin_m2 : np.ndarray
        Second partial derivative of the spin-s alms with respect to phi with
        factor 1/(sin theta)**2, shape (2, 12 * nside ** 2)
        if ret_slide is None, otherwise (2, len(ret_slice)).

    Notes
    -----
    The output maps are in the healpix format, with two entries: the real part
    and the imaginary part. The spin raising and lowering operators are applied
    to the spherical harmonic coefficients before the spin-s transform is
    computed. Note that the computation is done on the entire healpix grid,
    even if ret_slice is provided.

    """
    assert spin >= 0, spin
    assert hp.Alm.getlmax(gclm[0].size) == hp.Alm.getlmax(gclm[1].size)
    lmax = hp.Alm.getlmax(gclm[0].size)
    zbounds = np.sort(np.array(zbounds))
    # shape (2, 12 * nside ** 2),
    # first entry = real part, second entry imaginary part.

    # Computing the maps from the input alms,
    # where _sd[0] + 1j*_sd[1] is the spin (+s) fiedl,
    # and _sd[0] - 1j*_sd[1] is the spin (-s) field,
    # as well as the first derivatives with respect to theta and phi.
    _sd, d_dth, d_dphi_sin0 = alm2map_spin_der1(
        gclm,
        nside,
        spin,
        zbounds=zbounds,
        ret_slice=ret_slice
    )

    # Retrieving the application of two consecutive spin lowering operators
    _gclm = [
            hp.almxfl(
                alms,
                get_alpha_lower(spin, lmax) * get_alpha_lower(spin-1, lmax))
            for alms in gclm
        ]
    if spin > 2:
        _sm2_d = np.array(
            hp.alm2map_spin(_gclm, nside, spin - 2, lmax)
        )
    elif spin - 2 < 0:
        _sm2_d = np.array(
            hp.alm2map_spin(_gclm, nside, np.abs(spin - 2), lmax)
        )
        _sm2_d[1] *= -1
    else:
        _sm2_d = -np.array([hp.alm2map(alms, nside) for alms in _gclm])

    # Retrieving the application of two consecutive spin raising operators
    _gclm = [
        hp.almxfl(
            alms,
            get_alpha_raise(spin, lmax) * get_alpha_raise(spin+1, lmax))
        for alms in gclm
    ]
    _sp2_d = np.array(hp.alm2map_spin(_gclm, nside, spin + 2, lmax))

    # Retrieving the application of one spin raising
    # and one spin lowering operator
    _gclm = [
        hp.almxfl(
            alms,
            get_alpha_lower(spin, lmax) * get_alpha_raise(spin-1, lmax))
        for alms in gclm
    ]
    if spin == 0:
        _spm_d = -np.array([hp.alm2map(alms, nside) for alms in _gclm])
    else:
        _spm_d = np.array(hp.alm2map_spin(_gclm, nside, spin, lmax))

    # Retrieving the application of one spin lowering
    # and one spin raising operator
    _gclm = [
        hp.almxfl(
            alms,
            get_alpha_raise(spin, lmax) * get_alpha_lower(spin+1, lmax))
        for alms in gclm
    ]
    if spin == 0:
        _smp_d = -np.array([hp.alm2map(alms, nside) for alms in _gclm])
    else:
        _smp_d = np.array(hp.alm2map_spin(_gclm, nside, spin, lmax))

    # Computing the double partial derivative with respect to theta
    d2_dth2 = 0.25 * (_sp2_d + _sm2_d + _smp_d + _spm_d)

    # Computing the first term for the partial derivatives with respect to phi
    # and theta with the factor 1/sin(theta)
    # d_phi_sin0_d_th = 0.25 * (_sp2_d - _sm2_d)
    correction_1j = np.array([-1, 1]).reshape(2, 1)
    d_phi_sin0_d_th = - 0.25 * np.roll(_sp2_d - _sm2_d,
                                       axis=0,
                                       shift=1)*correction_1j
    # The roll is needed to account for the 1j factor and the correction_1j is
    # needed to account for the 1j * 1j contribution to the real part

    # Computing the first term for the double partial derivative with respect
    # to phi with the factor 1/sin(theta)**2
    d2_dphi2_sin_m2 = 0.25 * (_smp_d + _spm_d - (_sp2_d + _sm2_d))

    for iring in range(4 * nside - 1):
        startpix, nphi, kphi0, cth, sth = get_healpix_ring_pixel_layout(
            nside,
            iring
        )
        if zbounds[0] <= cth <= zbounds[1]:
            slic = slice(startpix, startpix + nphi)

            d_phi_sin0_d_th[:, slic] += np.roll(
                (spin * (cth / sth) ** 2 + spin/2.) * _sd
                - spin * (cth / sth) * d_dth,
                axis=0,
                shift=1
                )[:, slic]*correction_1j + ((cth / sth) * d_dphi_sin0)[:, slic]
            # The roll is needed to account for the 1j factor and
            # the correction_1j is needed to account for the 1j * 1j
            # contribution to the real part

            d2_dphi2_sin_m2[:, slic] -= (
                - spin**2 * (cth / sth) ** 2 * _sd
                + (cth / sth) * d_dth +
                np.roll(2 * spin * (cth / sth) * d_dphi_sin0, axis=0, shift=1)
                * correction_1j)[:, slic]
            # The roll is needed to account for the 1j factor
            #  and the correction_1j is needed to account for the 1j * 1j
            # contribution to the real part

    if ret_slice is not None:
        return (_sd[:, ret_slice],
                d_dth[:, ret_slice],
                d_dphi_sin0[:, ret_slice],
                d2_dth2[:, ret_slice],
                d_phi_sin0_d_th[:, ret_slice],
                d2_dphi2_sin_m2[:, ret_slice]
                )

    return _sd, d_dth, d_dphi_sin0, d2_dth2, d_phi_sin0_d_th, d2_dphi2_sin_m2
