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
    """Response coefficient of spin-s spherical harmonic to spin raising operator.

    Author: Julien Carron (j.carron@sussex.ac.uk)
    """
    ret = np.zeros(lmax + 1, dtype=float)
    ret[abs(s):] = np.sqrt(
        np.arange(abs(s) - s, lmax - s + 1) * np.arange(abs(s) + s + 1, lmax + s + 2)
    )
    return ret


def get_alpha_lower(s, lmax):
    """Response coefficient of spin-s spherical harmonic to spin lowering operator.

    Author: Julien Carron (j.carron@sussex.ac.uk)
    """
    ret = np.zeros(lmax + 1, dtype=float)
    ret[abs(s):] = -np.sqrt(
        np.arange(s + abs(s), lmax + s + 1) * np.arange(abs(s) - s + 1, lmax - s + 2)
    )
    return ret


def alm2map_spin_der1(gclm, nside, spin, zbounds=(-1.0, 1.0), ret_slice=None):
    """Returns spin-s transform '_{s}d' of alm,
    together with d/dtheta _{s}d and 1/sin tht d/dphi _{s}d.

    This crude version has three calls to spin-weight harmonics alm2map_spin.

    Author: Julien Carron (j.carron@sussex.ac.uk)
    """
    assert spin > 0, spin
    assert hp.Alm.getlmax(gclm[0].size) == hp.Alm.getlmax(gclm[1].size)
    lmax = hp.Alm.getlmax(gclm[0].size)
    zbounds = np.sort(np.array(zbounds))
    # shape (2, 12 * nside ** 2),
    # first entry = real part, second entry imaginary part.
    _sd = np.array(hp.alm2map_spin(gclm, nside, spin, lmax))
    if spin > 1:
        _gclm = [
            hp.almxfl(gclm[0], get_alpha_lower(spin, lmax)),
            hp.almxfl(gclm[1], get_alpha_lower(spin, lmax)),
        ]
        _sm1d = np.array(hp.alm2map_spin(_gclm, nside, spin - 1, lmax))
    else:
        _sm1d = -np.array(
            [
                hp.alm2map(
                    hp.almxfl(gclm[0], get_alpha_lower(spin, lmax)),
                    nside,
                    verbose=False,
                ),
                hp.alm2map(
                    hp.almxfl(gclm[1], get_alpha_lower(spin, lmax)),
                    nside,
                    verbose=False,
                ),
            ]
        )

    _gclm = [
        hp.almxfl(gclm[0], get_alpha_raise(spin, lmax)),
        hp.almxfl(gclm[1], get_alpha_raise(spin, lmax)),
    ]
    _sp1d = np.array(hp.alm2map_spin(_gclm, nside, spin + 1, lmax))

    d_dth = -0.5 * (_sp1d + _sm1d)

    d_dphi_sin0 = 0.5 * np.array([-_sp1d[1] + _sm1d[1], _sp1d[0] - _sm1d[0]])
    for iring in range(4 * nside - 1):
        startpix, nphi, kphi0, cth, sth = get_healpix_ring_pixel_layout(nside, iring)
        if zbounds[0] <= cth <= zbounds[1]:
            slic = slice(startpix, startpix + nphi)
            d_dphi_sin0[1, slic] -= spin * (cth / sth) * _sd[0, slic]
            d_dphi_sin0[0, slic] += spin * (cth / sth) * _sd[1, slic]
    if ret_slice is not None:
        return _sd[:, ret_slice], d_dth[:, ret_slice], d_dphi_sin0[:, ret_slice]

    return _sd, d_dth, d_dphi_sin0
