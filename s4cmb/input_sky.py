#!/usr/bin/python
"""
Script to simulate and handle input sky maps to be scanned.
Default file format is .fits containing healpix maps, and it comes with a
class HealpixFitsMap to handle it easily.
If you have a different I/O in your pipeline, just add a new class.

Author: Julien Peloton, j.peloton@sussex.ac.uk
"""
from __future__ import division, absolute_import, print_function

import glob
import os

import healpy as hp
import numpy as np
from astropy.io import fits as pyfits

class HealpixFitsMap():
    """ Class to handle fits file containing healpix maps """
    def __init__(self, input_filename,
                 do_pol=True, verbose=False,
                 no_ileak=False, no_quleak=False, ext_map_gal=False):
        """

        Parameters
        ----------
        input_filename : string
            Name of the file containing the sky maps.
        do_pol : bool, optional
            If True, load temperature and polarisation. Temperature only
            otherwise. Default is True.
        verbose : bool, optional
            If True, print out plenty of useless messages.
        no_ileak : bool, optional
            If True, load temperature and polarisation, but set the temperature
            to zero to avoid leakages.
        no_quleak : bool, optional
            If True, load temperature and polarisation, but set the
            polarisation to zero to avoid leakages.
        ext_map_gal : bool, optional
            Set it to True if you are reading a map in Galactic coordinate.
            (Planck maps for example).

        """
        self.input_filename = input_filename
        self.do_pol = do_pol
        self.verbose = verbose
        self.no_ileak = no_ileak
        self.no_quleak = no_quleak
        self.ext_map_gal = ext_map_gal

        self.nside = None
        self.I = None
        self.Q = None
        self.U = None

        self.load_healpix_fits_map()
        self.set_leakage_to_zero()

    def load_healpix_fits_map(self, force=False):
        """
        Load from disk into memory a sky map.

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
        >>> hpmap.load_healpix_fits_map(force=True)
        >>> print(hpmap.nside)
        16

        If the data is already loaded, it won't reload it by default
        >>> hpmap.load_healpix_fits_map()
        External data already present in memory
        """
        if self.nside is None or force:
            if self.do_pol:
                self.I, self.Q, self.U = hp.read_map(
                    self.input_filename, (0, 1, 2), verbose=self.verbose)
            else:
                self.I = hp.read_map(
                    self.input_filename, field=0, verbose=self.verbose)
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
        """
        ## Set temperature to zero to avoid I->QU leakage
        if self.no_ileak:
            self.I[:] = 0.0

        ## Set polarisation to zero to avoid QU leakage
        if self.no_quleak:
            if self.Q is not None:
                self.Q[:] = 0.0
            if self.U is not None:
                self.U[:] = 0.0

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
            lis[i] = ('HIERARCH ' + item[0], item[1], item[2])
        else:
            lis[i] = ('HIERARCH ' + item[0], item[1])
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
    theta_min = np.pi / 2. - ymax
    theta_max = np.pi / 2. - ymin
    fpix, lpix = hp.ang2pix(nside, [theta_min, theta_max], [0., 2.*np.pi])
    pixs = np.arange(fpix, lpix + 1, dtype=np.int)

    theta, phi = hp.pix2ang(nside, pixs)
    if xmin < 0:
        phi[phi > np.pi] = (phi[phi > np.pi] - 2 * np.pi)
    good = (theta >= theta_min) * (theta <= theta_max) * \
        (phi <= xmax) * (phi >= xmin)
    obspix = pixs[good]
    obspix.sort()

    return obspix

def create_sky_map(cl_fn, nside=16, seed=548397):
    """
    Create full sky map from input cl.

    Parameters
    ----------
    cl_fn : string
        Name of the file containing cl (CAMB lensed cl format)
    nside : int, optional
        Resolution for the output map.

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
    ell, TT, EE, BB, TE = np.loadtxt(cl_fn).T

    ## Take out the normalisation...
    llp = ell * (ell + 1.) / (2 * np.pi)

    np.random.seed(seed)
    I, Q, U = hp.synfast([TT / llp, EE / llp, BB / llp, TE / llp], nside,
                         lmax=2*nside, mmax=None, alm=False,
                         pol=True, pixwin=False,
                         fwhm=0.0, sigma=None, new=True,
                         verbose=False)
    return I, Q, U

def write_healpix_cmbmap(output_filename, data, nside, fits_IDL=False,
                         coord=None, colnames=['I', 'Q', 'U'], nest=False):
    """
    Write healpix fits map in full sky mode or custom partial sky,
    i.e. file with obspix and CMB_fields. Input data have to be a list
    with n fields to be written.

    / ! \
    By default, even the full sky mode write the maps in partial mode, in
    the sense that the data is compressed. so unless you know what you
    are doing, always choose partial_custom=False.
    / ! \

    Parameters
    ----------
    output_filename : string
        Name of the output file (.fits).
    data : list of 1d array(s)
        Data to save on disk.
    nside : int
        Resolution of the map. Must be a power of 2.
    fits_IDL : bool
        If True, store the data reshaped in row of 1024 (IDL style).
        Default is False.
    coord : string
        The system of coordinates in which the data are
        (G(alactic), C(elestial), and so on). Default is None.
    colnames : list of strings
        The name of each data vector to be saved.
    nest : bool, optional
        If True, save the data in the nest scheme. Default is False (i.e.
        data are saved in the RING format).

    Examples
    ----------
    >>> nside = 16
    >>> I, Q, U = np.random.rand(3, hp.nside2npix(nside))
    >>> colnames = ['I', 'Q', 'U']
    >>> write_healpix_cmbmap('myfits_to_test_.fits',
    ...     data=[I, Q, U], nside=nside, colnames=colnames)
    """
    ## Write the header
    extra_header = []
    for c in colnames:
        extra_header.append(('column_names', c))
    extra_header = add_hierarch(extra_header)

    hp.write_map(output_filename, data, fits_IDL=fits_IDL,
                 coord=coord, column_names=None, partial=True,
                 extra_header=extra_header)

def write_dummy_map(filename='myfits_to_test_.fits', nside=16):
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
    colnames = ['I', 'Q', 'U']
    write_healpix_cmbmap(filename, data=[I, Q, U],
                         nside=nside, colnames=colnames)

def remove_test_data(has_id='_to_test_', silent=True):
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
    fns = glob.glob('*' + has_id + '*')
    if not silent:
        print('Removing files: ', fns)
    for fn in fns:
        os.remove(fn)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
    remove_test_data(has_id='_to_test_', silent=True)
