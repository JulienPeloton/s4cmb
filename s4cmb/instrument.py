#!/usr/bin/python
"""
Script to simulate the hardware of a CMB experiment.
* focal plane
* pointing model parameters of the telescope
* beam parameters of the bolometers
* polarisation angle of the bolometers

Author: Julien Peloton, j.peloton@sussex.ac.uk
"""
from __future__ import division, absolute_import, print_function

import os
import copy
import glob
import datetime
import numpy as np

def coordinates_on_grid(pix_size=None, row_size=None,
                        nx=None, nx2=None,
                        max_points=None):
    """
    Return the x and y coordinates of points on a grid.
    The grid is centered on (0, 0).

    Parameters
    ----------
    pix_size : float, optional
        Size of each pixel. User should either provide
        `pix_size` or `row_size` (but not both at the same time).
    row_size : float, optional
        Size of each row. User should either provide
        `pix_size` or `row_size` (but not both at the same time).
    nx : int, optional
        Number of pixels per row/column. User should either provide
        `nx` or `nx2` (but not both at the same time).
    nx2 : int, optional
        Total number of pixels in the array. User should either provide
        `nx` or `nx2` (but not both at the same time).
    max_points : int, optional
        If nx2 is specified, `max_points` defines the maximum number of points
        to return. If None, set to `nx2`.

    Returns
    ----------
    coordinates : ndarray (2, nx[:max_points] * nx[:max_points])
        x and y coordinates of the pixels.

    Examples
    ----------
    Make a grid with 2 points per row, spaced by 1 unit
    >>> coordinates_on_grid(pix_size=1., nx=2)
    ... # doctest: +NORMALIZE_WHITESPACE
    array([[-0.5,  0.5, -0.5,  0.5],
           [-0.5, -0.5,  0.5,  0.5]])

    Same grid, but specifying the total number of points
    >>> coordinates_on_grid(pix_size=1., nx2=4)
    ... # doctest: +NORMALIZE_WHITESPACE
    array([[-0.5,  0.5, -0.5,  0.5],
           [-0.5, -0.5,  0.5,  0.5]])

    You can also specify the maximum number of points to return
    >>> coordinates_on_grid(pix_size=1., nx2=4, max_points=3)
    ... # doctest: +NORMALIZE_WHITESPACE
    array([[-0.5,  0.5, -0.5],
           [-0.5, -0.5,  0.5]])

    If you specify a total number of points which does not fit totally inside
    a squared grid, then it will return only the point you ask on a bigger grid
    Note that it works as if you would have specify max_points for a higher nx2
    >>> coordinates_on_grid(pix_size=1., nx2=3)
    ... # doctest: +NORMALIZE_WHITESPACE
    array([[-0.5,  0.5, -0.5],
           [-0.5, -0.5,  0.5]])

    You should specify either the number of pixel per row column (nx),
    or the total number of point in the grid but not both at the same time:
    >>> coordinates_on_grid(pix_size=1., nx=2, nx2=4)
    ... # doctest: +NORMALIZE_WHITESPACE, +ELLIPSIS
    Traceback (most recent call last):
      ...
    AssertionError: You should specify either the number of pixel
    per row column (nx), or the total number of point in the grid (nx2).

    Idem for pix_size and row_size
    >>> coordinates_on_grid(pix_size=1., row_size=4., nx=2)
    ... # doctest: +NORMALIZE_WHITESPACE, +ELLIPSIS
    Traceback (most recent call last):
      ...
    AssertionError: You should specify either the size of
    a pixel (pix_size), or the size of a row (row_size).
    """
    if (nx is None and nx2 is None) or (nx is not None and nx2 is not None):
        raise AssertionError('You should specify either the ' +
                             'number of pixel per row column (nx), or ' +
                             'the total number of point in the grid (nx2).\n')

    if (pix_size is None and row_size is None) or \
            (pix_size is not None and row_size is not None):
        raise AssertionError('You should specify either the ' +
                             'size of a pixel (pix_size), ' +
                             'or the size of a row (row_size).\n')

    if nx2 is not None:
        ## Look for the closest number with square root being an integer
        nx2_tmp = copy.copy(nx2)
        while True:
            nx = np.sqrt(nx2_tmp)
            if int(nx) == nx:
                nx = int(nx)
                break
            else:
                nx2_tmp += 1
    else:
        nx2 = nx**2

    if max_points is None:
        max_points = nx2

    if pix_size is None:
        pix_size = row_size / nx
    elif row_size is None:
        row_size = pix_size * nx

    ix1 = np.arange(nx)
    xs1 = (ix1 - (nx - 1.) / 2.) * pix_size
    x2, y2 = np.meshgrid(xs1, xs1)

    coordinates = np.array(
        (x2.flatten()[:max_points],
         y2.flatten()[:max_points]))

    return coordinates

def convert_pair_to_bolometer_position(xcoord_pairs, ycoord_pairs):
    """
    Return the position of bolometers given the position of pairs.

    Parameters
    ----------
    xcoord_pairs : ndarray
        Array of length `npair` containing the coordinate
        of the pairs of bolometers along the x axis.
    ycoord_pairs : ndarray
        Array of length `npair` containing the coordinate
        of the pairs of bolometers along the y axis.

    Returns
    ----------
    xcoord_bolometers : ndarray
        Array of length `2 * npair` containing the coordinate of
        the bolometers along the x axis.
    ycoord_bolometers : ndarray
        Array of length `2 * npair` containing the coordinate of
        the bolometers along the y axis.

    Examples
    ----------
    >>> fp = FocalPlane(verbose=False)
    >>> xp, yp = coordinates_on_grid(row_size=fp.fp_size, nx2=4)
    >>> print(xp, yp)
    [-15.  15. -15.  15.] [-15. -15.  15.  15.]
    >>> xb, yb = convert_pair_to_bolometer_position(xp, yp)
    >>> print(xb, yb) # doctest: +NORMALIZE_WHITESPACE
    [-15. -15.  15.  15. -15. -15.  15.  15.]
    [-15. -15. -15. -15.  15.  15.  15.  15.]
    """
    nbolometer = 2 * len(xcoord_pairs)
    xcoord_bolometers = np.dstack(
        (xcoord_pairs, xcoord_pairs)).reshape((1, nbolometer))[0]
    ycoord_bolometers = np.dstack(
        (ycoord_pairs, ycoord_pairs)).reshape((1, nbolometer))[0]

    return xcoord_bolometers, ycoord_bolometers

def show_focal_plane(bolo_xcoord, bolo_ycoord, bolo_polangle=None,
                     fn_out='plot_hardware_map_test.png',
                     save_on_disk=True, display=False):
    """
    Show the focal plane of the instrument, split in two panels:
    top and bottom bolometers

    Parameters
    ----------
    bolo_xcoord : 1d array
        Bolometers x coordinates in the focal plane.
    bolo_ycoord : 1d array
        Bolometers y coordinates in the focal plane.
    bolo_polangle : 1d array, optional
        Bolometer intrinsic polarisation angle orientation. If provided,
        it is used to color code the figure.
    fn_out : string, optional
        Name of the output file containing the plot of the focal plane.
        Provide the extension (format: png or pdf).
    save_on_disk : bool
        If True, save the plot on disk.
    display : bool
        If True, show the plot.

    Examples
    ---------
    >>> fp = FocalPlane(verbose=False)
    >>> show_focal_plane(fp.bolo_xcoord, fp.bolo_ycoord, fp.bolo_polangle,
    ...     save_on_disk=False, display=False)
    """
    if not display:
        import matplotlib as mpl
        mpl.use('Agg')
        import matplotlib.pyplot as pl
        pl.ioff()
    else:
        import matplotlib.pyplot as pl

    if bolo_polangle is None:
        bolo_polangle = np.ones_like(bolo_xcoord)

    fig, ax = pl.subplots(1, 2, figsize=(10, 5))
    ## Top pixel
    ax[0].scatter(bolo_xcoord[::2], bolo_ycoord[::2],
                  c=bolo_polangle[::2], alpha=1, s=30, cmap=pl.cm.jet)
    ax[0].scatter(bolo_xcoord[::2], bolo_ycoord[::2],
                  c='black', s=30, marker='|',
                  label='Top pixel', alpha=0.6)
    ax[0].set_ylabel('y position (cm)')
    ax[0].set_xlabel('x position (cm)')
    ax[0].set_title('Top pixels')

    ## Bottom pixel
    ax[1].scatter(bolo_xcoord[1::2], bolo_ycoord[1::2],
                  c=bolo_polangle[1::2], alpha=1, s=30, cmap=pl.cm.jet)
    ax[1].scatter(bolo_xcoord[1::2], bolo_ycoord[1::2],
                  c='black', s=30, marker='_',
                  label='Bottom pixel', alpha=0.6)
    ax[1].set_ylabel('y position (cm)')
    ax[1].set_xlabel('x position (cm)')
    ax[1].set_title('Bottom pixels')

    if save_on_disk:
        pl.savefig(fn_out)
    if display:
        pl.show()
    pl.clf()

def convert_cm_to_rad(xcm, ycm, conversion):
    """
    Convert positions in cm of the pairs of bolometers or bolometers
    in the focal plane into positions in radian.

    Parameters
    ----------
    xcm : ndarray
        x coordinates of pairs of bolometers or bolometers in cm.
    ycm : ndarray
        y coordinates of pairs of bolometers or bolometers in cm.
    conversion : float
        Conversion factor in rad/cm.

    Returns
    ----------
    xrad : ndarray
        x coordinates of pairs of bolometers or bolometers in radians.
        It has the same length as `xcm`.
    yrad : ndarray
        y coordinates of pairs of bolometers or bolometers in radians.
        It has the same length as `ycm`.

    Examples
    ----------
    Focal plane of 60 cm diameter and mirror
    giving a 3 deg projection on the sky by default
    >>> fp = FocalPlane(verbose=False)
    >>> projected_fp_size = 3. ## degrees
    >>> xcm, ycm = coordinates_on_grid(
    ...                 row_size=fp.fp_size, nx2=fp.npair)
    >>> print(convert_cm_to_rad(xcm, ycm,
    ...     conversion=projected_fp_size / fp.fp_size * np.pi / 180.))
    ...     # doctest: +NORMALIZE_WHITESPACE
    (array([-0.01308997,  0.01308997, -0.01308997,  0.01308997]),
     array([-0.01308997, -0.01308997,  0.01308997,  0.01308997]))

    """
    return np.array(xcm) * conversion, np.array(ycm) * conversion

def construct_beammap(beamprm, ct, cb, nx, pix_size):
    """
    Construct the pixel beam maps
    (sum and difference of bolometer beam maps)

    Parameters
    ----------
    beamprm : beam_model instance
        Instance of beam_model.
    ct : int
        Index of the top bolometer in the pair.
    cb : int
        Index of the bottom bolometer in the pair.
    nx : int
        Number of pixels per row/column (in pixel).
    pix_size : float
        Size of each pixel (in radian).

    Returns
    ----------
    summap : ndarray
        Beam map made of the sum of bolometer beam maps.
    diffmap : ndarray
        Beam map made of the difference of bolometer beam maps.

    Examples
    ----------
    Note that bolometers within the same pixel are neighbour bolometers
    that is (ct, cb) = (0, 1) for example.
    >>> fp = FocalPlane(verbose=False)
    >>> bm = BeamModel(fp, verbose=False)
    >>> pix_size = 0.5 / 180. * np.pi / 60. # 0.5 arcmin in rad
    >>> summap, diffmap = construct_beammap(
    ...     beamprm=bm, ct=0, cb=1, nx=4, pix_size=pix_size)
    >>> print(summap) # doctest: +NORMALIZE_WHITESPACE
    [[ 0.77520676  0.86809111  0.86809111  0.77520676]
     [ 0.86809111  0.97210474  0.97210474  0.86809111]
     [ 0.86809111  0.97210474  0.97210474  0.86809111]
     [ 0.77520676  0.86809111  0.86809111  0.77520676]]
    """
    # Translate beams to origin and maintain differential pointing
    dx = beamprm.xpos[ct] - beamprm.xpos[cb]
    dy = beamprm.ypos[ct] - beamprm.ypos[cb]

    tx = 0.5 * dx
    bx = -0.5 * dx
    ty = 0.5 * dy
    by = -0.5 * dy

    xy2f = coordinates_on_grid(pix_size=pix_size, nx=nx)

    tmap = gauss2d(xy2f, tx, ty,
                   beamprm.Amp[ct], beamprm.sig_1[ct],
                   beamprm.sig_2[ct],
                   beamprm.ellip_ang[ct]).reshape((nx, nx))

    bmap = gauss2d(xy2f, bx, by,
                   beamprm.Amp[cb], beamprm.sig_1[cb],
                   beamprm.sig_2[cb],
                   beamprm.ellip_ang[cb]).reshape((nx, nx))

    summap = 0.5 * (tmap + bmap)
    diffmap = 0.5 * (tmap - bmap)

    return summap, diffmap

def gauss2d(xy, x_0, y_0, Amp, sig_xp, sig_yp, psi):
    """
    2D Gaussian model for beams.

    Parameters
    ----------
    xy : 2d array
        Columns are projected coordinates (xproj,yproj)
    x_0: float
        x coordinate of the center of the Gaussian.
    y_0: float
        y coordinate of the center of the Gaussian.
    Amp: float
        Amplitude of the Gaussian.
    sig_xp: float
        Sigma for the Gaussian in x' coordinate system (rotated
        system if psi != 0).
    sig_yp: float
        Sigma for the Gaussian in y' coordinate system (rotated
        system if psi != 0).
    psi: float
        Angle between normal coordinate system and
        primed system (normal is primed if psi = 0).

    Returns
    ----------
    z : 1d array
        Flatten beam map of shape (xy[1].shape, )

    Examples
    ----------
    >>> fp = FocalPlane(verbose=False)
    >>> bm = BeamModel(fp, verbose=False)
    >>> pix_size = 0.5 / 180. * np.pi / 60. # 0.5 arcmin in rad
    >>> xy = coordinates_on_grid(pix_size=pix_size, nx=4)
    >>> gauss2d(xy, x_0=0, y_0=0, Amp=1.,
    ...     sig_xp=bm.sig_1[0],
    ...     sig_yp=bm.sig_2[0], psi=0)
    ... # doctest: +NORMALIZE_WHITESPACE
    array([ 0.77520676,  0.86809111,  0.86809111,  0.77520676,
            0.86809111,  0.97210474,  0.97210474,  0.86809111,
            0.86809111,  0.97210474,  0.97210474,  0.86809111,
            0.77520676,  0.86809111,  0.86809111,  0.77520676])
    """

    x_1 = xy[0, :] - x_0
    y_1 = xy[1, :] - y_0

    xy_1 = np.array([x_1, y_1])

    psi2 = -psi * np.pi / 180.0

    # x/y coordinates make an angle psi with the xp/yp input coordinates
    R = np.array(
        [[np.cos(psi2), -np.sin(psi2)], [np.sin(psi2), np.cos(psi2)]])
    p = np.dot(R, xy_1)

    u = p[0, :]**2 / (2 * sig_xp**2) + p[1, :]**2 / (2 * sig_yp**2)

    # Hide underflow by clipping beam function at -430dB level
    mask = u < 100
    z = Amp * np.exp(-u * mask) * mask

    return z

class Hardware():
    """ Class to load all the hardware and models of the instrument in once """
    def __init__(self,
                 ncrate=1, ndfmux_per_crate=1, nsquid_per_mux=1,
                 npair_per_squid=4, fp_size=60.,
                 fwhm=3.5, beam_seed=58347,
                 projected_fp_size=3.,
                 pm_name='5params',
                 type_hwp='CRHWP', freq_hwp=2., angle_hwp=0., verbose=False):
        """
        This class creates the data used to model the instrument:
        * focal plane
        * pointing model parameters of the telescope
        * beam parameters of the bolometers
        * polarisation angle of the bolometers

        Parameters
        ----------
        ncrate : int
            Number of crate plate.
        ndfmux_per_crate : int
            Number of MUX board per crate.
        nsquid_per_mux : int
            Number of SQUID per MUX.
        npair_per_squid : int
            Number of pair of bolometers per SQUID.
        fp_size : float, optional
            The size of the focal plane in cm. Default is 60 cm.
        fwhm : float, optional
            Full Width Half Maximum of the beam (in arcmin).
            Default = 3.5 arcmin.
        beam_seed : int
            Seed used to generate angle of rotation of beam axes.
            Default is 58347.
        projected_fp_size : float, optional
            Size of the focal plane on the sky (in degree). This has to
            do with the size of the mirror. Default = 3 degrees.
        pm_name : string, optional
            The pointing model to load. Currently, only the five-parameter
            pointing model (Mangum 2001) is implemented (pm_name = 5params).
        verbose : boolean, optional
            If True, print out a number of useful comments for verboseging.

        Examples
        ----------
        >>> instrument = Hardware()
        """
        self.focal_plane = FocalPlane(ncrate, ndfmux_per_crate,
                                      nsquid_per_mux, npair_per_squid,
                                      fp_size, verbose)

        self.beam_model = BeamModel(self.focal_plane, fwhm, beam_seed,
                                    projected_fp_size, verbose)

        self.pointing_model = PointingModel(pm_name)

        self.half_wave_plate = HalfWavePlate(type_hwp, freq_hwp, angle_hwp)

class FocalPlane():
    """ Class to handle the focal plane of the instrument. """
    def __init__(self,
                 ncrate=1, ndfmux_per_crate=1, nsquid_per_mux=1,
                 npair_per_squid=4, fp_size=60., verbose=False):
        """
        Initialise our focal plane.

        Note.
        The polarisation angle model consists in defining detector polarisation
        angle. The focal plane is cut in quadrants (Crate).
        Within a quadrant, pixels are categorized into two: Q and U pixels.
        Q and U pixels have 45 degrees difference in their polarisation angle,
        and form lines within quadrant. Each pixel contains a top and a
        bottom bolometer, with 90 degrees difference in the polarisation angle.
        Then, you go from one quadrant to another by a global 90 deg rotation
        of the polarisation angle starting with quadrant 0 having
        theta_{Q, top} = 0 deg.

        Parameters
        ----------
        ncrate : int
            Number of crate plate.
        ndfmux_per_crate : int
            Number of MUX board per crate.
        nsquid_per_mux : int
            Number of SQUID per MUX.
        npair_per_squid : int
            Number of pair of bolometers per SQUID.
        fp_size : float, optional
            The size of the focal plane in cm. Default is 60 cm.
        verbose : boolean, optional
            If True, print out a number of useful comments for verboseging.
        """
        self.ncrate = ncrate
        self.ndfmux_per_crate = ndfmux_per_crate
        self.nsquid_per_mux = nsquid_per_mux
        self.npair_per_squid = npair_per_squid

        ## Total number of pairs and bolometers in the focal plane
        self.npair = self.ncrate * self.ndfmux_per_crate * \
            self.nsquid_per_mux * self.npair_per_squid
        self.nbolometer = self.npair * 2

        self.fp_size = fp_size

        self.verbose = verbose

        self.make_focal_plane()

    def make_focal_plane(self):
        """
        Create the hardware map of the instrument,
        taht is the focal plane geometry,
        the bolometers id, the wiring, and so on.
        The terminology used here is taken from the Polarbear experiment.
        The hierarchy is the following:

        +-------------------------------------------------------------------+
        |CRATE -> DFMUX -> SQUID -> BOLOMETER (top & bottom)
        |  |        |        |          |
        |  v        v        v          v
        |  id       id       id         id, xCoordinate, yCoordinate,
        |                               IndexInFocalPlane, polangle_orientation
        |                               IndexInSquid.
        +--------------------------------------------------------------------+

        Examples
        ----------
        >>> fp = FocalPlane(verbose=True)
        Hardware map generated...
        """
        ## Retrieve coordinate of the pairs inside the focal plane
        # xcoord, ycoord = self.compute_pairs_coordinates(self.npair)
        xcoord, ycoord = coordinates_on_grid(row_size=self.fp_size,
                                             nx2=self.npair)

        ## Initialise
        self.crate_id, self.dfmux_id = [], []
        self.squid_id, self.bolo_id = [], []
        self.bolo_index_in_squid, self.bolo_index_in_fp = [], []
        self.bolo_xcoord, self.bolo_ycoord, self.bolo_polangle = [], [], []
        ## Construct the hardware map
        max_hit = False
        while max_hit is False:
            bolo_index = 0
            pair_index = 0
            squid_index = 0
            dfmux_index = 0
            for crate in range(self.ncrate):
                ## CRATE
                self.crate_id.append('Cr{:03}'.format(crate))

                for dfmux in range(self.ndfmux_per_crate):
                    ## DFMUX
                    self.dfmux_id.append('Cr{:03}Df{:03}'.format(
                            crate, dfmux_index))

                    for squid in range(self.nsquid_per_mux):
                        ## SQUID
                        self.squid_id.append('Cr{:03}Df{:03}Sq{:03}'.format(
                            crate, dfmux_index, squid_index))

                        for pair in range(self.npair_per_squid):
                            ## BOLOMETER

                            ## Split Q/U
                            boloQ = 2*pair
                            boloU = 2*pair + 1

                            if int(pair_index/np.sqrt(self.npair)) % 2 == 0:
                                shift = 0.
                                shiftinv = 45.
                            else:
                                shift = 45.
                                shiftinv = 0.

                            ## Top pixel
                            ## Position of the bolometer within the SQUID
                            self.bolo_index_in_squid.append(boloQ)
                            self.bolo_index_in_fp.append(bolo_index)
                            self.bolo_id.append(
                                'Cr{:03}Df{:03}Sq{:03}Bo{:03}t'.format(
                                    crate, dfmux_index, squid_index, boloQ))
                            self.bolo_xcoord.append(xcoord[pair_index])
                            self.bolo_ycoord.append(ycoord[pair_index])

                            ## Q/U pixels
                            # bolo in pairs are separated by 90 deg,
                            # and each quadran is separated by 90 deg.
                            if bolo_index % 4 == 0:
                                shift_ = 45.
                                shiftinv_ = 0.
                            else:
                                shift_ = 0.
                                shiftinv_ = 45.
                            if xcoord[pair_index] < 0 and \
                                    ycoord[pair_index] >= 0:
                                angle = shift
                            elif xcoord[pair_index] >= 0 and \
                                    ycoord[pair_index] >= 0:
                                angle = 90. + shift_
                            elif xcoord[pair_index] >= 0 and \
                                    ycoord[pair_index] <= 0:
                                angle = 180. + shiftinv
                            elif xcoord[pair_index] <= 0 and \
                                    ycoord[pair_index] <= 0:
                                angle = 270. + shiftinv_

                            self.bolo_polangle.append(angle)

                            ## Move in bolo space (2 bolos/pair)
                            bolo_index += 1

                            ## Bottom pixel
                            ## Top pixel
                            ## Position of the bolometer within the SQUID
                            self.bolo_index_in_squid.append(boloU)
                            self.bolo_index_in_fp.append(bolo_index)
                            self.bolo_id.append(
                                'Cr{:03}Df{:03}Sq{:03}Bo{:03}b'.format(
                                    crate, dfmux_index, squid_index, boloU))
                            self.bolo_xcoord.append(xcoord[pair_index])
                            self.bolo_ycoord.append(ycoord[pair_index])

                            ## 90 degree difference wrt top bolometer
                            self.bolo_polangle.append((angle + 90) % 360)

                            ## Move again in bolo space (2 bolos/pair)
                            bolo_index += 1

                            ## Move in pair space
                            pair_index += 1

                            ## Close the job if you hit the maximum number of
                            ## bolometers or pairs.
                            try:
                                assert bolo_index < self.nbolometer, \
                                    'Hardware map generated...'
                                assert pair_index < self.npair, \
                                    'Hardware map generated...'
                            except AssertionError as e:
                                if self.verbose:
                                    print(str(e))
                                max_hit = True
                                break
                        squid_index += 1
                    dfmux_index += 1

    def get_indices(self, name='Cr'):
        """
        Returns Cr(ate), Df(mux) or Sq(uid) indices.

        Parameters
        ----------
        name : string
            Which component you want the indices.
            Should be in ['Cr', 'Sq', 'Df'].
            Note that if you want bolometer indices, just use
            self.bolo_index_in_fp (for global indices in the focal plane) or
            self.bolo_index_in_squid (for local indices in SQUIDs).

        Returns
        ----------
        indices : list of int
            List containing the indices of the component for all bolometers.

        Examples
        ----------
        >>> fp = FocalPlane()
        >>> squid_indices = fp.get_indices('Sq')
        >>> print(squid_indices)
        [0, 0, 0, 0, 0, 0, 0, 0]

        """
        assert name in ['Cr', 'Sq', 'Df'], \
            ValueError("name must be in ['Cr', 'Sq', 'Df'].")

        indices = [
            int(comp.split(name)[-1][:3]) for comp in self.bolo_id]

        return indices

class BeamModel():
    """ Class to handle the beams of the detectors """
    def __init__(self,
                 focal_plane, fwhm=3.5, beam_seed=58347,
                 projected_fp_size=3., verbose=False):
        """
        Parameters
        ----------
        focal_plane : focal_plane instance
            Instance of focal_plane containing focal plane parameters.
        fwhm : float, optional
            Full Width Half Maximum of the beam (in arcmin).
            Default = 3.5 arcmin.
        beam_seed : int
            Seed used to generate angle of rotation of beam axes.
            Default is 58347.
        projected_fp_size : float, optional
            Size of the focal plane on the sky (in degree). This has to
            do with the size of the mirror. Default = 3 degrees.
        verbose : boolean, optional
            If True, print out a number of useful comments for verboseging.
        """
        ## Focal plane parameters
        self.focal_plane = focal_plane

        ## Beam model and mirror parameters
        self.fwhm = fwhm
        self.beam_seed = beam_seed
        self.projected_fp_size = projected_fp_size

        ## Paths and names
        self.verbose = verbose

        self.beamprm = self.generate_beam_parameters()

    def generate_beam_parameters(self):
        """
        Construct the beam parameters such as position of the centroids,
        ellipticity, angle of rotation, and associated errors.

        Returns
        ----------
        beamprm_fields : dictionary
            Dictionary containing beam parameters

        Examples
        ----------
        >>> fp = FocalPlane(verbose=False)
        >>> bm = BeamModel(fp, verbose=False)
        >>> bm.generate_beam_parameters()
        >>> print(bm.xpos) # doctest: +NORMALIZE_WHITESPACE
        [-0.01308997 -0.01308997  0.01308997  0.01308997
         -0.01308997 -0.01308997  0.01308997  0.01308997]
        """

        beamprm_header = ['Amp', 'Amp_err',
                          'ellip_ang', 'ellip_ang_err',
                          'sig_1', 'sig_1_err',
                          'sig_2', 'sig_2_err',
                          'xpos', 'xpos_err',
                          'ypos', 'ypos_err']

        for b in beamprm_header:
            setattr(self, b, np.zeros(self.focal_plane.nbolometer))

        ## Position of the bolometers
        ## (bolometers cm -> bolometers radians)
        self.xpos, self.ypos = convert_cm_to_rad(
            self.focal_plane.bolo_xcoord, self.focal_plane.bolo_ycoord,
            conversion=self.projected_fp_size /
            self.focal_plane.fp_size * np.pi / 180.)

        ## Generate Gaussian beams.
        ## FWHM arcmin -> FWHM rad -> sigma rad
        FWHM_rad = self.fwhm / 60. * np.pi / 180.
        sigma_rad = FWHM_rad / np.sqrt(8 * np.log(2))

        self.sig_1 = np.ones(self.focal_plane.nbolometer) * sigma_rad
        self.sig_2 = np.ones(self.focal_plane.nbolometer) * sigma_rad

        ## Angle of rotation for the ellipses.
        ## Between -90 and 90 degrees.
        state = np.random.RandomState(self.beam_seed)
        self.ellip_ang = state.uniform(-90, 90, self.focal_plane.nbolometer)

        ## Amplitude of the beams
        ## Default is one.
        self.Amp = np.ones(self.focal_plane.nbolometer)

class PointingModel():
    """ Class to handle the pointing model of the telescope """
    def __init__(self, pm_name='5params'):
        """
        We focus on a five-parameter pointing model (Mangum 2001) to
        characterize the relationship between the telescope's encoder
        readings and its true boresight pointing on the sky.
        The parameters described in this reference are

        * IA, the azimuth encoder zero offset,
        * IE, the elevation encoder zero offset,
        * CA, the collimation error of the electromagnetic axis,
        * AN, the azimuth axis offset/misalignment (north-south) and
        * AW, the azimuth offset/misalignment (east-west).

        Parameters
        ----------
        pm_name : string, optional
            The pointing model to load. Currently, only the five-parameter
            pointing model (Mangum 2001) is implemented (pm_name = 5params).

        Examples
        --------
        >>> pm = PointingModel()
        >>> [(i, round(j,3)) for i, j in zip(
        ...     pm.allowed_params.split(), pm.value_params)]
        ... # doctest: +NORMALIZE_WHITESPACE
        [('ia', -10.285), ('ie', 8.74),
         ('ca', -15.598), ('an', -0.51),
         ('aw', 0.109)]

        >>> pm = PointingModel(pm_name='super-model-trop-bien')
        ... # doctest: +NORMALIZE_WHITESPACE, +ELLIPSIS
        Traceback (most recent call last):
         ...
        ValueError: Only the five-parameter pointing model
        (Mangum 2001) is implemented for the moment (pm_name = 5params)
        """
        self.pm_name = pm_name

        if self.pm_name == '5params':
            self.five_parameter_pointing_model()
        else:
            raise ValueError('Only the five-parameter ' +
                             'pointing model (Mangum 2001) is implemented ' +
                             'for the moment (pm_name = 5params)')

    def five_parameter_pointing_model(self):
        """
        Parameters based on Polarbear configuration.
        """
        self.allowed_params = 'ia ie ca an aw'

        self.value_params = np.array([-10.28473073, 8.73953334,
                                      -15.59771781, -0.50977716, 0.10858016])

        ## Set this to zero for the moment
        self.RMS_AZ = 0.0
        self.RMS_EL = 0.0
        self.RMS_RESID = 0.0

class HalfWavePlate():
    """ Class to handle the Half-Wave Plate (HWP) """
    def __init__(self, type_hwp='CRHWP', freq_hwp=2., angle_hwp=0.):
        """
        This class provides routines to compute the HWP angles.
        This can be use later to generate timestreams.

        Parameters
        ----------
        type_hwp : string, optional
            The type of HWP that you want to mount on your instrument.
            * CRWHP: continously rotating HWP.
            * stepped: stepped HWP (once a CES).
        freq_hwp : float, optional
            The frequency of rotation of the HWP in Hz.
        angle_hwp : float, optional
            The offset of the HWP in degree.

        """
        self.type_hwp = type_hwp
        self.freq_hwp = freq_hwp
        self.angle_hwp = angle_hwp

        if self.type_hwp not in ['CRHWP', 'stepped']:
            raise ValueError("`type_hwp` has to be 'CRHWP' or 'stepped'.")

        if self.type_hwp is 'stepped' and freq_hwp != 0.0:
            raise AssertionError("You cannot have a stepped HWP and non-" +
                                 "zero frequency! set freq_hwp=0.0 " +
                                 "if you want a stepped HWP.")

    def compute_HWP_angles(self, sample_rate=1., size=1):
        """
        Generate HWP angles which can be use later to generate timestreams.

        Parameters
        ----------
        sample_rate : float, optional
            Sample rate of the detectors
        size : int, optional
            Length of the vector of angles (number of time samples).

        Returns
        ----------
        HWP_angles : ndarray
            1d array of size `size` containting the HWP angles in radian.

        Examples
        Continously rotating HWP at 2 Hz starting at 0 degree.
        >>> hwp = HalfWavePlate(type_hwp='CRHWP', freq_hwp=2., angle_hwp=0.)
        >>> hwp.compute_HWP_angles(sample_rate=100., size=10)
        ... # doctest: +NORMALIZE_WHITESPACE
        array([ 0.        ,  0.12566371,  0.25132741,  0.37699112,  0.50265482,
            0.62831853,  0.75398224,  0.87964594,  1.00530965,  1.13097336])

        Stepped HWP with 30 degrees angle
        >>> hwp = HalfWavePlate(type_hwp='stepped',
        ...     freq_hwp=0.0, angle_hwp=30.)
        >>> hwp.compute_HWP_angles(sample_rate=100., size=10)
        ... # doctest: +NORMALIZE_WHITESPACE
        array([ 0.52359878,  0.52359878,  0.52359878,  0.52359878,  0.52359878,
            0.52359878,  0.52359878,  0.52359878,  0.52359878,  0.52359878])

        For a stepped HWP, the frequency must be zero
        >>> hwp = HalfWavePlate(type_hwp='stepped',
        ...     freq_hwp=1.0, angle_hwp=30.)
        ... # doctest: +NORMALIZE_WHITESPACE, +ELLIPSIS
        Traceback (most recent call last):
            ...
        AssertionError: You cannot have a stepped HWP and non-zero frequency!
        set freq_hwp=0.0 if you want a stepped HWP.
        """
        angle = self.angle_hwp * np.pi / 180.

        HWP_angles = np.array(
            [angle + t * (self.freq_hwp / sample_rate) *
             2. * np.pi for t in range(size)])

        return HWP_angles

    def update_hardware(self, new_type_hwp, new_freq_hwp, new_angle_hwp):
        """
        Change the behaviour of the HWP.

        Parameters
        ----------
        new_type_hwp : string, optional
            The type of HWP that you want to mount on your instrument.
            * CRWHP: continously rotating HWP.
            * stepped: stepped HWP (once a CES).
        new_freq_hwp : float, optional
            The frequency of rotation of the HWP in Hz.
        new_angle_hwp : float, optional
            The offset of the HWP in degree.

        Examples
        ----------
        Continously rotating HWP at 2 Hz starting at 0 degree.
        >>> hwp = HalfWavePlate(type_hwp='CRHWP', freq_hwp=2., angle_hwp=0.)
        >>> hwp.compute_HWP_angles(sample_rate=100., size=10)
        ... # doctest: +NORMALIZE_WHITESPACE
        array([ 0.        ,  0.12566371,  0.25132741,  0.37699112,  0.50265482,
            0.62831853,  0.75398224,  0.87964594,  1.00530965,  1.13097336])

        For some reason, our HWP died!
        >>> hwp.update_hardware(new_type_hwp='stepped',
        ...     new_freq_hwp=0.0, new_angle_hwp=0.0)
        >>> hwp.compute_HWP_angles(sample_rate=100., size=10)
        ... # doctest: +NORMALIZE_WHITESPACE
        array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.])

        """
        self.type_hwp = new_freq_hwp
        self.freq_hwp = new_freq_hwp
        self.angle_hwp = new_angle_hwp


if __name__ == "__main__":
    import doctest
    doctest.testmod()
