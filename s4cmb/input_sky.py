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
Script to simulate and handle input sky maps to be scanned.
Default file format is .fits containing healpix maps, and it comes with a
class HealpixFitsMap to handle it easily.
If you have a different I/O in your pipeline, just add a new class.

Author: Julien Peloton, peloton@lal.in2p3.fr
        Giulio Fabbian, g.fabbian@sussex.ac.uk
"""
from __future__ import division, absolute_import, print_function

import glob
import os
import warnings

import healpy as hp
import numpy as np

from s4cmb.tools import alm2map_spin_der1


class HealpixFitsMap:
    """ Class to handle fits file containing healpix maps """

    def __init__(
        self,
        input_filename,
        do_pol=True,
        verbose=False,
        fwhm_in=0.0,
        fwhm_in2=None,
        nside_in=16,
        lmax=None,
        map_seed=53543,
        no_ileak=False,
        no_quleak=False,
        compute_derivatives=None,
        derivatives_type="T",
        ext_map_gal=False,
    ):
        """

        Parameters
        ----------
        input_filename : string, or list of strings
            Either fits file containing the sky maps (data will just be
            loaded), or CAMB lensed cl file (.dat) containing lensed
            power spectra with order ell, TT, EE, BB, TE (maps will be
            created on-the-fly), or a list of 3 fits files containing alms
            (maps will be created on-the-fly).
        do_pol : bool, optional
            If True, load temperature and polarisation. Temperature only
            otherwise. Default is True.
        verbose : bool, optional
            If True, print out plenty of useless messages.
        fwhm_in : float, optional
            If input_filename is a CAMB lensed cl file, the generated maps will
            be convolved with a beam having this fwhm_in. In arcmin.
            No effect if you provide maps directly.
        fwhm_in2 : float, optional
            If provided, will generate another set of I, Q, U with this
            resolution (useful for dichroic detectors). Default is None.
            No effect if you provide maps directly.
        nside_in : int, optional
            If input_filename is a CAMB lensed cl file, the maps will be
            generated at a resolution nside_in. No effect if you provide maps
            directly.
        lmax : None or int, optional
            Maximum multipole when creating a map from cl. If none, it
            is set automatically to 2*nside_in
        map_seed : int, optional
            If input_filename is a CAMB lensed cl file, this is the seed used
            to create the maps. No effect if you provide maps directly.
        no_ileak : bool, optional
            If True, load temperature and polarisation, but set the temperature
            to zero to avoid leakages.
        no_quleak : bool, optional
            If True, load temperature and polarisation, but set the
            polarisation to zero to avoid leakages.
        ext_map_gal : bool, optional
            Set it to True if you are reading a map in Galactic coordinate.
            (Planck maps for example).
        compute_derivatives : bool, optional
            If True, return derivatives of relevant Stokes component according
            to the derivatives_type option.
            Be sure that you have enough memory!
        derivatives_type : str, optional
            If 'T' ('P') present in the string return 1st and 2nd derivatives
            of the input temperature (Q and U) map (t=theta, p=phi).
            Note that d/dp is already divided by sin(theta).
            Options:
            - 'T1' will compute 1st derivatives dI/dt, dI/dp alone
            - 'P1' will compute 1st derivatives dQ/dt, dU/dt, dQ/dp, dU/dp.
            Be sure that you have enough memory!
        """
        self.input_filename = input_filename
        self.do_pol = do_pol
        self.verbose = verbose
        self.no_ileak = no_ileak
        self.no_quleak = no_quleak
        self.ext_map_gal = ext_map_gal
        self.fwhm_in = fwhm_in
        self.fwhm_in2 = fwhm_in2
        self.nside_in = nside_in
        if lmax is None:
            self.lmax = 2 * self.nside_in
        else:
            self.lmax = lmax
        self.map_seed = map_seed
        self.compute_derivatives = compute_derivatives
        self.derivatives_type = derivatives_type
        self.I = None
        self.Q = None
        self.U = None

        self.I2 = None
        self.Q2 = None
        self.U2 = None

        fromalms = False
        if type(self.input_filename) == list:
            if self.verbose:
                print("Reading sky maps from alms file...")
            self.load_healpix_fits_map_from_alms()
            fromalms = True
        elif self.input_filename[-4:] == ".dat":
            if self.verbose:
                print("Creating sky maps from cl file...")
            self.create_healpix_fits_map()
        elif self.input_filename[-5:] == ".fits":
            if self.verbose:
                print("Reading sky maps from fits file...")
            self.load_healpix_fits_map()
        else:
            raise IOError(
                """Input file not understood! Should be either a
                fits file containing the sky maps
                (data will be loaded), or a CAMB lensed cl file
                (.dat) containing lensed power spectra with
                order ell, TT, EE, BB, TE
                (maps will be created on-the-fly)."""
            )

        self.set_leakage_to_zero()

        if self.compute_derivatives:
            if "T" in self.derivatives_type:
                self.compute_intensity_derivatives(fromalm=fromalms)
            if "P" in self.derivatives_type:
                self.compute_pol_derivatives(fromalm=fromalms)

    def load_healpix_fits_map(self, force=False):
        """
        Load from disk into memory a sky map.

        Not updated for dichroic for the moment.

        Parameters
        ----------
        force : bool
            If true, force to load the maps in memory even if it is already
            loaded. Default is False.

        Examples
        ----------
        Let's generate fake data
        >>> filename = 'myfits_to_test_.fits'
        >>> write_dummy_map(filename, nside=16)

        Let's now read the data
        >>> hpmap = HealpixFitsMap(input_filename=filename)
        >>> print(hpmap.nside)
        16

        If the data is already loaded, it won't reload it by default
        >>> hpmap.load_healpix_fits_map()
        External data already present in memory

        But you can force it
        >>> hpmap.load_healpix_fits_map(force=True)
        """
        if self.I is None or force:
            if self.do_pol:
                self.I, self.Q, self.U = hp.read_map(
                    self.input_filename, (0, 1, 2), verbose=self.verbose
                )
            else:
                self.I = hp.read_map(
                    self.input_filename, field=0, verbose=self.verbose
                )
            self.nside = hp.npix2nside(len(self.I))
        else:
            print("External data already present in memory")

    def load_healpix_fits_map_from_alms(self, force=False):
        """
        Load from disk into memory alms and make sky maps.

        Parameters
        ----------
        force : bool
            If true, force to load the maps in memory even if it is already
            loaded. Default is False.

        Examples
        ----------
        Let's generate fake data
        >>> np.random.seed(548397)
        >>> sky_maps = create_sky_map('s4cmb/data/test_data_set_lensedCls.dat')
        >>> alms = hp.map2alm(sky_maps)
        >>> filenames = ['myalms_to_test_tlm.fits', 'myalms_to_test_elm.fits',
        ...     'myalms_to_test_blm.fits']
        >>> for fn, alm in zip(filenames, alms):
        ...     hp.write_alm(fn, alm)

        Let's now read the data
        >>> hpmap = HealpixFitsMap(input_filename=filenames)
        >>> print(hpmap.nside)
        16

        If the data is already loaded, it won't reload it by default
        >>> hpmap.load_healpix_fits_map_from_alms()
        External data already present in memory

        But you can force it
        >>> hpmap.load_healpix_fits_map_from_alms(force=True)

        You can also generate 2 sets of maps with different resolution which
        is useful for dichroic detectors
        >>> hpmap = HealpixFitsMap(input_filename=filenames, fwhm_in=3.5,
        ...     fwhm_in2=1.8, nside_in=16,)
        >>> hasattr(hpmap, 'Q2')
        True
        """
        if self.I is None or force:
            if self.do_pol:
                tlm = hp.read_alm(self.input_filename[0])
                elm = hp.read_alm(self.input_filename[1])
                blm = hp.read_alm(self.input_filename[2])

                self.I, self.Q, self.U = hp.alm2map(
                    [tlm, elm, blm],
                    nside=self.nside_in,
                    pixwin=False,
                    fwhm=self.fwhm_in / 60.0 * np.pi / 180.0,
                    sigma=None,
                    pol=True,
                    inplace=False,
                    verbose=self.verbose,
                )
                if self.fwhm_in2 is not None:
                    self.I2, self.Q2, self.U2 = hp.alm2map(
                        [tlm, elm, blm],
                        nside=self.nside_in,
                        pixwin=False,
                        fwhm=self.fwhm_in2 / 60.0 * np.pi / 180.0,
                        sigma=None,
                        pol=True,
                        inplace=False,
                        verbose=self.verbose,
                    )
            else:
                tlm = hp.read_alm(self.input_filename[0])

                self.I = hp.alm2map(
                    tlm,
                    nside=self.nside_in,
                    pixwin=False,
                    fwhm=self.fwhm_in / 60.0 * np.pi / 180.0,
                    sigma=None,
                    pol=False,
                    inplace=False,
                    verbose=self.verbose,
                )
                if self.fwhm_in2 is not None:
                    self.I2 = hp.alm2map(
                        tlm,
                        nside=self.nside_in,
                        pixwin=False,
                        fwhm=self.fwhm_in2 / 60.0 * np.pi / 180.0,
                        sigma=None,
                        pol=False,
                        inplace=False,
                        verbose=self.verbose,
                    )

            self.nside = hp.npix2nside(len(self.I))
        else:
            print("External data already present in memory")

    def create_healpix_fits_map(self, force=False):
        """
        Create sky maps from cl file.
        Do nothing if data already presents in the memory.

        Parameters
        ----------
        force : bool
            If true, force to recreate the maps in memory even
            if it is already loaded. Default is False.

        Examples
        ----------
        Let's generate the map from a CAMB file
        >>> filename = 's4cmb/data/test_data_set_lensedCls.dat'
        >>> hpmap = HealpixFitsMap(input_filename=filename, fwhm_in=3.5,
        ...     nside_in=16, map_seed=489237)
        >>> print(hpmap.nside)
        16

        If the data is already loaded, it won't reload it by default
        >>> hpmap.create_healpix_fits_map()
        External data already present in memory

        But you can force it
        >>> hpmap.create_healpix_fits_map(force=True)

        You can also load 2 sets of maps with different resolution, which
        is useful for dichroic detectors
        >>> filename = 's4cmb/data/test_data_set_lensedCls.dat'
        >>> hpmap = HealpixFitsMap(input_filename=filename, fwhm_in=3.5,
        ...     fwhm_in2=1.8, nside_in=16, map_seed=489237)
        >>> hasattr(hpmap, 'I2')
        True
        """
        if self.I is None or force:
            if self.do_pol:
                self.I, self.Q, self.U = create_sky_map(
                    self.input_filename,
                    nside=self.nside_in,
                    FWHM=self.fwhm_in,
                    seed=self.map_seed,
                    lmax=self.lmax,
                )
                if self.fwhm_in2 is not None:
                    self.I2, self.Q2, self.U2 = create_sky_map(
                        self.input_filename,
                        nside=self.nside_in,
                        FWHM=self.fwhm_in2,
                        seed=self.map_seed,
                        lmax=self.lmax,
                    )
            else:
                self.I = create_sky_map(
                    self.input_filename,
                    nside=self.nside_in,
                    FWHM=self.fwhm_in,
                    seed=self.map_seed,
                    lmax=self.lmax,
                )
                if self.fwhm_in2 is not None:
                    self.I2 = create_sky_map(
                        self.input_filename,
                        nside=self.nside_in,
                        FWHM=self.fwhm_in2,
                        seed=self.map_seed,
                        lmax=self.lmax,
                    )
            self.nside = hp.npix2nside(len(self.I))
        else:
            print("External data already present in memory")

    def set_leakage_to_zero(self):
        """
        Remove either I, Q or U to remove possible leakages

        Examples
        ----------
        Test with no input intensity
        >>> write_dummy_map('myfits_to_test_.fits')
        >>> hpmap = HealpixFitsMap('myfits_to_test_.fits', no_ileak=True)
        >>> print(hpmap.I)
        [ 0.  0.  0. ...,  0.  0.  0.]

        Test with no input polarisation
        >>> write_dummy_map('myfits_to_test_.fits')
        >>> hpmap = HealpixFitsMap('myfits_to_test_.fits', no_quleak=True)
        >>> print(hpmap.Q, hpmap.U)
        [ 0.  0.  0. ...,  0.  0.  0.] [ 0.  0.  0. ...,  0.  0.  0.]

        If you have two sets of maps, it will remove leakages from the two sets
        >>> filename = 's4cmb/data/test_data_set_lensedCls.dat'
        >>> hpmap = HealpixFitsMap(input_filename=filename, fwhm_in=3.5,
        ...     fwhm_in2=1.8, nside_in=16, map_seed=489237,
        ...     no_ileak=True, no_quleak=True)
        >>> print(hpmap.I, hpmap.I2)
        [ 0.  0.  0. ...,  0.  0.  0.] [ 0.  0.  0. ...,  0.  0.  0.]
        """
        # Set temperature to zero to avoid I->QU leakage
        if self.no_ileak:
            self.I[:] = 0.0
            if self.I2 is not None:
                self.I2[:] = 0.0

        # Set polarisation to zero to avoid QU leakage
        if self.no_quleak:
            if self.Q is not None:
                self.Q[:] = 0.0
            if self.U is not None:
                self.U[:] = 0.0

            if self.Q2 is not None:
                self.Q2[:] = 0.0
            if self.U2 is not None:
                self.U2[:] = 0.0

    def compute_intensity_derivatives(self, fromalm=False):
        """
        Compute derivatives of the input temperature map (healpix).
        Not updated for dichroic for the moment.

        Parameters
        ----------
        fromalm : bool, optional
            If True, loads alm file from disk instead of fourier
            transform the input map. Automatically turns True if you input
            alm files. False otherwise.

        Examples
        ----------
        >>> filename = 's4cmb/data/test_data_set_lensedCls.dat'
        >>> hpmap = HealpixFitsMap(input_filename=filename, fwhm_in=3.5,
        ...     nside_in=16, compute_derivatives=True, map_seed=489237)
        >>> hasattr(hpmap, 'dIdp')
        True

        """
        if fromalm:
            alm = hp.read_alm(self.input_filename[0])
        else:
            alm = hp.map2alm(self.I, self.lmax)

        # lmax = hp.Alm.getlmax(alm.size)
        if "T1" in self.derivatives_type:
            junk, self.dIdt, self.dIdp = hp.alm2map_der1(
                alm, self.nside_in, self.lmax
            )
        else:
            # computes first and second derivative as derivatives of spin-1
            # transform of a scalar field with _1Elm=sqrt(l(l+1))Ilm _1Blm=0
            l = np.arange(self.lmax + 1)
            grad = np.sqrt(l * (l + 1))
            curl = np.zeros_like(alm)
            dervs = alm2map_spin_der1([hp.almxfl(alm, grad), curl], self.nside_in, 1)
            self.dIdt = dervs[0][0]
            self.dIdp = dervs[0][1]
            self.d2Id2t = dervs[1][0]
            self.d2Id2p = dervs[2][1]
            self.d2Idpdt = dervs[2][0]

    def compute_pol_derivatives(self, fromalm=False):
        """
        Compute derivatives of the input polarization components (healpix).
        Not updated for dichroic for the moment.

        Parameters
        ----------
        fromalm : bool, optional
            If True, loads alm file from disk instead of fourier
            transform the input map. Automatically turns True if you input
            alm files. False otherwise.

        Examples
        ----------
        >>> filename = 's4cmb/data/test_data_set_lensedCls.dat'
        >>> hpmap = HealpixFitsMap(input_filename=filename, fwhm_in=3.5,
        ...     nside_in=16, compute_derivatives=True, map_seed=489237)
        >>> hasattr(hpmap, 'dIdp')
        True

        """
        if fromalm:
            Elm = hp.read_alm(self.input_filename[1])
            Blm = hp.read_alm(self.input_filename[2])
        else:
            alm = hp.map2alm([self.I, self.Q, self.U], self.lmax)
            Elm = alm[1]
            Blm = alm[2]
        # lmax = hp.Alm.getlmax(Elm.size)
        if "P1" in self.derivatives_type:
            out = alm2map_spin_der1([Elm, Blm], self.nside_in, 2)
            self.dQdt = out[1][0]
            self.dUdt = out[1][1]
            self.dQdp = out[2][0]
            self.dUdp = out[2][1]
        else:
            warnings.warn("""
            Computation of second order polarization derivatives not
            implemented yet. Set to 0.
            """)

            out = alm2map_spin_der1([Elm, Blm], self.nside_in, 2)
            self.dQdt = out[1][0]
            self.dUdt = out[1][1]
            self.dQdp = out[2][0]
            self.dUdp = out[2][1]
            self.d2Qd2t = np.zeros_like(self.dQdt)
            self.d2Qd2p = np.zeros_like(self.dQdt)
            self.d2Qdpdt = np.zeros_like(self.dQdt)
            self.d2Ud2t = np.zeros_like(self.dQdt)
            self.d2Ud2p = np.zeros_like(self.dQdt)
            self.d2Udpdt = np.zeros_like(self.dQdt)


def add_hierarch(lis):
    """
    Convert in correct format for fits header.

    Parameters
    ----------
    lis: list of tuples
        Contains tuples (keyword, value [, comment]).

    Returns
    ----------
    lis : list of strings
        Contains strings in the pyfits header format.

    Examples
    ----------
    >>> lis = [['toto', 3, 'I am a comment']]
    >>> add_hierarch(lis)
    [('HIERARCH toto', 3, 'I am a comment')]
    """
    for i, item in enumerate(lis):
        if len(item) == 3:
            lis[i] = ("HIERARCH " + item[0], item[1], item[2])
        else:
            lis[i] = ("HIERARCH " + item[0], item[1])
    return lis


def get_obspix(xmin, xmax, ymin, ymax, nside):
    """
    Given RA/Dec boundaries, return the observed pixels in the healpix scheme.

    Parameters
    ----------
    xmin : float
        x coordinate of the bottom left corner in radian. (RA min)
    xmax : float
        x coordinate of the bottom right corner in radian. (RA max)
    ymin : float
        y coordinate of the top left corner in radian. (Dec min)
    ymax : float
        y coordinate of the top right corner in radian. (Dec max)
    nside : int
        Resolution of the healpix map.

    Returns
    ----------
    obspix : 1d array of int
        The list of observed pixels.

    Examples
    ----------
    >>> get_obspix(-np.pi/2, np.pi/2,
    ...     -np.pi/2, np.pi/2, nside=2) # doctest: +NORMALIZE_WHITESPACE
    array([ 0,  3,  4,  5, 10, 11, 12, 13, 14,
           18, 19, 20, 21, 26, 27, 28, 29,
           30, 34, 35, 36, 37, 42, 43, 44])
    """
    theta_min = np.pi / 2.0 - ymax
    theta_max = np.pi / 2.0 - ymin
    fpix, lpix = hp.ang2pix(nside, [theta_min, theta_max], [0.0, 2.0 * np.pi])
    pixs = np.arange(fpix, lpix + 1, dtype=np.int)

    theta, phi = hp.pix2ang(nside, pixs)
    if xmin < 0:
        phi[phi > np.pi] = phi[phi > np.pi] - 2 * np.pi
    good = (theta >= theta_min) * (theta <= theta_max) * (phi <= xmax) * (phi >= xmin)
    obspix = pixs[good]
    obspix.sort()

    return obspix


def LamCyl(ra, dec):
    """
    Referred to as cylindrical equal-area in the USGS report, assuming
    that the parallel of true scale is zero
    """
    return ra, np.sin(dec)


def SFL(ra, dec):
    """SFL stands for Sanson-Flamsteed. In the USGS report this is
    referred to as the Sinusoidal Projection. It is equal-area. Parallels
    are equally spaced straight lines. Scale is true along central meridian
    and all paralles."""
    return ra * np.cos(dec), dec


def deSFL(x, y):
    return x / np.cos(y), y


def deLamCyl(x, y):
    return x, np.arcsin(y)


def create_sky_map(cl_fn, nside=16, FWHM=0.0, seed=548397, lmax=None):
    """
    Create full sky map from input cl.

    Parameters
    ----------
    cl_fn : string
        Name of the file containing cl (CAMB lensed cl format)
    nside : int, optional
        Resolution for the output map.
    FWHM : float
        The fwhm of the Gaussian used to smooth the map (applied on alm).
        In arcmin.

    Returns
    ----------
    maps : ndarray
        Maps of the sky (I, Q, U) of size 12 * nside**2.

    Examples
    ----------
    Create a sky map. Seed is fixed for testing purposes.
    >>> np.random.seed(548397)
    >>> sky_maps = create_sky_map('s4cmb/data/test_data_set_lensedCls.dat')
    >>> print(sky_maps[0])
    [  55.51567033   50.94330727   39.69851524 ...,   36.2265932   107.64964085
       80.8613084 ]
    """
    if lmax is None:
        lmax = 2 * nside
    ell, TT, EE, BB, TE = np.loadtxt(cl_fn).T

    # Take out the normalisation...
    llp = ell * (ell + 1.0) / (2 * np.pi)

    # Arcmin to rad
    FWHM_rad = FWHM / 60.0 * np.pi / 180.0

    np.random.seed(seed)
    I, Q, U = hp.synfast(
        [TT / llp, EE / llp, BB / llp, TE / llp],
        nside,
        lmax=lmax,
        mmax=None,
        alm=False,
        pol=True,
        pixwin=False,
        fwhm=FWHM_rad,
        sigma=None,
        new=True,
        verbose=False,
    )
    return I, Q, U


def write_healpix_cmbmap(
    output_filename,
    data,
    fits_IDL=False,
    coord=None,
    colnames=["I", "Q", "U"],
    partial=True,
    nest=False,
):
    """
    Write healpix fits map in full sky mode or partial sky,
    Input data have to be a list with n fields to be written.

    Parameters
    ----------
    output_filename : string
        Name of the output file (.fits).
    data : list of 1d array(s)
        Data to save on disk.
    fits_IDL : bool
        If True, store the data reshaped in row of 1024 (IDL style).
        Default is False.
    coord : string
        The system of coordinates in which the data are
        (G(alactic), C(elestial), and so on). Default is None.
    colnames : list of strings
        The name of each data vector to be saved.
    partial : bool
        If True, store only non-zero pixels. Default is True.
    nest : bool, optional
        If True, save the data in the nest scheme. Default is False (i.e.
        data are saved in the RING format).

    Examples
    ----------
    >>> nside = 16
    >>> I, Q, U = np.random.rand(3, hp.nside2npix(nside))
    >>> colnames = ['I', 'Q', 'U']
    >>> write_healpix_cmbmap('myfits_to_test_.fits',
    ...     data=[I, Q, U], colnames=colnames)
    """
    # Write the header
    extra_header = []
    for c in colnames:
        extra_header.append(("column_names", c))
    extra_header = add_hierarch(extra_header)

    hp.write_map(
        output_filename,
        data,
        fits_IDL=fits_IDL,
        coord=coord,
        column_names=None,
        partial=partial,
        extra_header=extra_header,
        overwrite=True,
    )


def write_dummy_map(filename="myfits_to_test_.fits", nside=16):
    """
    Write dummy file on disk for test purposes.

    Parameters
    ----------
    filename : string, optional
        Name of the output file (.fits)
    nside : int
        Resolution of the maps.

    Examples
    ----------
    >>> write_dummy_map()
    """
    nside = 16
    I, Q, U = np.random.rand(3, hp.nside2npix(nside))
    colnames = ["I", "Q", "U"]
    write_healpix_cmbmap(filename, data=[I, Q, U], colnames=colnames)


def remove_test_data(has_id="_to_test_", silent=True):
    """
    Remove data with name containing the `has_id`.

    Parameters
    ----------
    has_id : string
        String included in filename(s) to remove.

    Examples
    ----------
    >>> file = open('file_to_erase_.txt', 'w')
    >>> file.close()
    >>> remove_test_data(has_id='_to_erase_', silent=False)
    Removing files:  ['file_to_erase_.txt']
    """
    fns = glob.glob("*" + has_id + "*")
    if not silent:
        print("Removing files: ", fns)
    for fn in fns:
        os.remove(fn)


if __name__ == "__main__":
    import doctest

    if np.__version__ >= "1.14.0":
        np.set_printoptions(legacy="1.13")
    doctest.testmod()
    remove_test_data(has_id="_to_test_", silent=True)
