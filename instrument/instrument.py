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

## I/O
from astropy.io import fits as pyfits
import xml.etree.ElementTree as ET

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

class hardware():
    """ Class to load the hardware and models of the instrument """
    def __init__(self,
                 ncrate=1, ndfmux_per_crate=1, nsquid_per_mux=1,
                 npair_per_squid=4, fp_size=60.,
                 FWHM=3.5, beam_seed=58347,
                 projected_fp_size=3.,
                 pm_name='5params',
                 type_HWP='CRHWP', freq_HWP=2., angle_HWP=0.,
                 output_folder='./', name='_to_test_', debug=False):
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
        FWHM : float, optional
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
        output_folder : string, optional
            The folder where the data will be stored.
        name : string, optional
            Tag for the output file (without the extension).
        debug : boolean, optional
            If True, print out a number of useful comments for debugging.

        Examples
        ----------
        >>> instrument = hardware()
        """
        self.focal_plane = focal_plane(ncrate, ndfmux_per_crate,
                                       nsquid_per_mux, npair_per_squid,
                                       fp_size, output_folder, name, debug)

        self.beam_model = beam_model(self.focal_plane, FWHM, beam_seed,
                                     projected_fp_size, output_folder,
                                     name, debug)

        self.pointing_model = pointing_model(pm_name, output_folder, name)

        self.polarisation_angle_model = polarisation_angle_model(
            self.focal_plane, output_folder, name)

        self.half_wave_plate = half_wave_plate(type_HWP, freq_HWP, angle_HWP)

class focal_plane():
    """ Class to handle the focal plane of the instrument. """
    def __init__(self,
                 ncrate=1, ndfmux_per_crate=1, nsquid_per_mux=1,
                 npair_per_squid=4, fp_size=60.,
                 output_folder='./', name='_to_test_', debug=False):
        """
        Initialise our focal plane and save the configuration
        inside a xml file.

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
        output_folder : string, optional
            The folder where the xml file will be stored.
        name : string, optional
            The name of the xml file (without the extension).
        debug : boolean, optional
            If True, print out a number of useful comments for debugging.
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
        self.output_folder = output_folder
        self.name = name
        self.output_file = os.path.join(
            self.output_folder, 'focal_plane_' + self.name + '.xml')

        self.debug = debug

        self.create_hwmap()

    def create_hwmap(self):
        """
        Create the hardware map of the instrument,
        that is the xml file containing the focal plane geometry,
        the bolometers id, the wiring, and so on.
        The terminology used here is taken from the Polarbear experiment.
        The hierarchy is the following:

        +-------------------------------------------------------------------+
        |CRATE -> DFMUX -> SQUID -> BOLOMETER (top & bottom)
        |  |        |        |          |
        |  v        v        v          v
        |  id       id       id         id, xCoordinate, yCoordinate,
        |                               focalPlaneIndex, polangle_orientation
        |                               polarizationMode,
        |                               polarizationOrientation, channel
        +--------------------------------------------------------------------+

        Examples
        ----------
        >>> fp = focal_plane(debug=True)
        Hardware map generated...
        Hardware map written at ./focal_plane__to_test_.xml
        """
        ## Retrieve coordinate of the pairs inside the focal plane
        # xcoord, ycoord = self.compute_pairs_coordinates(self.npair)
        xcoord, ycoord = coordinates_on_grid(row_size=self.fp_size,
                                             nx2=self.npair)

        HWmap = ET.Element('HardwareMap')

        ## Construct the hardware map
        max_hit = False
        while max_hit is False:
            bolo_id = 0
            pair_id = 0
            polangle = []
            name_bolo = []
            for crate in range(self.ncrate):
                ## CRATE
                ET.SubElement(HWmap, 'Crate')
                HWmap[crate].set('id', 'Cr%d' % crate)

                for dfmux in range(self.ndfmux_per_crate):
                    ## DFMUX
                    ET.SubElement(HWmap[crate], 'DfMuxBoard')
                    HWmap[crate][dfmux].set('id', 'Cr%dDf%d' % (crate, dfmux))

                    for squid in range(self.nsquid_per_mux):
                        ## SQUID
                        ET.SubElement(HWmap[crate][dfmux], 'Squid')
                        HWmap[crate][dfmux][squid].set('id', 'Cr%dDf%dSq%d' % (
                            crate, dfmux, squid))

                        for pair in range(self.npair_per_squid):
                            ## BOLOMETER

                            ## Split Q/U
                            boloQ = 2*pair
                            boloU = 2*pair + 1

                            if int(pair_id/np.sqrt(self.npair)) % 2 == 0:
                                shift = 0.
                                shiftinv = 45.
                            else:
                                shift = 45.
                                shiftinv = 0.

                            ## Top pixel
                            ET.SubElement(
                                HWmap[crate][dfmux][squid], 'Bolometer')
                            HWmap[crate][dfmux][squid][boloQ].set(
                                'focalPlaneIndex', '1')
                            HWmap[crate][dfmux][squid][boloQ].set(
                                'polarizationMode', '1')
                            ## Position of the bolometer within the SQUID
                            HWmap[crate][dfmux][squid][boloQ].set(
                                'channel', '%d' % boloQ)

                            name_bolo_tmp = 'Cr%dDf%dSq%d_%dt' % (
                                crate, dfmux, squid, pair_id)
                            name_bolo.append(name_bolo_tmp)
                            HWmap[crate][dfmux][squid][boloQ].set(
                                'polarizationOrientation', '0')
                            HWmap[crate][dfmux][squid][boloQ].set(
                                'id', name_bolo_tmp)
                            HWmap[crate][dfmux][squid][boloQ].set(
                                'xCoordinate', '%.4f' % xcoord[pair_id])
                            HWmap[crate][dfmux][squid][boloQ].set(
                                'yCoordinate', '%.4f' % ycoord[pair_id])

                            ## Q/U pixels
                            # bolo in pairs are separated by 90 deg,
                            # and each quadran is separated by 90 deg.
                            if bolo_id % 4 == 0:
                                shift_ = 45.
                                shiftinv_ = 0.
                            else:
                                shift_ = 0.
                                shiftinv_ = 45.
                            if xcoord[pair_id] < 0 and ycoord[pair_id] >= 0:
                                angle = shift
                            elif xcoord[pair_id] >= 0 and ycoord[pair_id] >= 0:
                                angle = 90. + shift_
                            elif xcoord[pair_id] >= 0 and ycoord[pair_id] <= 0:
                                angle = 180. + shiftinv
                            elif xcoord[pair_id] <= 0 and ycoord[pair_id] <= 0:
                                angle = 270. + shiftinv_
                            HWmap[crate][dfmux][squid][boloQ].set(
                                'polangle_orientation', '%.2f' % (angle))
                            polangle.append(angle)

                            ## Move in bolo space (2 bolos/pair)
                            bolo_id += 1

                            ## Bottom pixel
                            ET.SubElement(
                                HWmap[crate][dfmux][squid], 'Bolometer')
                            HWmap[crate][dfmux][squid][boloU].set(
                                'focalPlaneIndex', '1')
                            HWmap[crate][dfmux][squid][boloU].set(
                                'polarizationMode', '1')
                            ## Position of the bolometer within the SQUID
                            HWmap[crate][dfmux][squid][boloU].set(
                                'channel', '%d' % boloU)

                            name_bolo_tmp = 'Cr%dDf%dSq%d_%db' % (
                                crate, dfmux, squid, pair_id)
                            name_bolo.append(name_bolo_tmp)
                            HWmap[crate][dfmux][squid][boloU].set(
                                'polarizationOrientation', '1')
                            HWmap[crate][dfmux][squid][boloU].set(
                                'id', name_bolo_tmp)
                            HWmap[crate][dfmux][squid][boloU].set(
                                'xCoordinate', '%.4f' % xcoord[pair_id])
                            HWmap[crate][dfmux][squid][boloU].set(
                                'yCoordinate', '%.4f' % ycoord[pair_id])

                            ## Q/U pixels
                            if xcoord[pair_id] < 0 and ycoord[pair_id] >= 0:
                                angle = 90.+shift
                            elif xcoord[pair_id] >= 0 and ycoord[pair_id] >= 0:
                                angle = 180.+shift_
                            elif xcoord[pair_id] >= 0 and ycoord[pair_id] <= 0:
                                angle = 270.+shiftinv
                            elif xcoord[pair_id] <= 0 and ycoord[pair_id] <= 0:
                                angle = shiftinv_
                            HWmap[crate][dfmux][squid][boloU].set(
                                'polangle_orientation', '%.2f' % (angle))
                            polangle.append(angle)
                            ## Move again in bolo space (2 bolos/pair)
                            bolo_id += 1

                            ## Move in pair space
                            pair_id += 1

                            ## Close the job if you hit the maximum number of
                            ## bolometers or pairs.
                            try:
                                assert bolo_id < self.nbolometer, \
                                    'Hardware map generated...'
                                assert pair_id < self.npair, \
                                    'Hardware map generated...'
                            except AssertionError as e:
                                if self.debug:
                                    print(str(e))
                                max_hit = True
                                break

        tree = ET.ElementTree(HWmap)
        tree.write(self.output_file)
        if self.debug:
            print('Hardware map written at {}'.format(self.output_file))

    @staticmethod
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
        >>> fp = focal_plane(debug=False)
        >>> xp, yp = coordinates_on_grid(row_size=fp.fp_size, nx2=4)
        >>> print(xp, yp)
        [-15.  15. -15.  15.] [-15. -15.  15.  15.]
        >>> xb, yb = fp.convert_pair_to_bolometer_position(xp, yp)
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

    @staticmethod
    def unpack_hwmap(fn, tag, key, dtype=str):
        """
        Routine to extract focal plane data from the hardware map generated
        by create_hwmap(). Particularly useful for interfacing.

        Parameters
        ----------
        fn : string
            xml file containing the hardware map written by create_hwmap.
        tag : string
            Tag corresponding to the layer in the xml file (Crate, DfMuxBoard,
            Squid, Bolometer).
        key : string
            The key corresponding to the data to be extracted.
            `key` should be in `tag`.
        dtype : type
            Type of the data to unpack. Default is string.

        Returns
        ----------
        data : ndarray
            The data corresponding to `key`.

        Examples
        ----------
        >>> fp = focal_plane(debug=False)

        Return the id of the Crate boards in the focal plane (one here)
        >>> fp.unpack_hwmap(fn='./focal_plane__to_test_.xml',
        ...     tag='Crate', key='id') # doctest: +NORMALIZE_WHITESPACE
        array(['Cr0'], dtype='|S3')

        Return the id of the DfMux boards in the focal plane (one here)
        >>> fp.unpack_hwmap(fn='./focal_plane__to_test_.xml',
        ...     tag='DfMuxBoard', key='id') # doctest: +NORMALIZE_WHITESPACE
        array(['Cr0Df0'], dtype='|S6')

        Return the id of the Squids in the focal plane (one here)
        >>> fp.unpack_hwmap(fn='./focal_plane__to_test_.xml',
        ...     tag='Squid', key='id') # doctest: +NORMALIZE_WHITESPACE
        array(['Cr0Df0Sq0'], dtype='|S9')

        Return the id of the 8 bolometers (4 pairs) in the focal plane
        >>> fp.unpack_hwmap(fn='./focal_plane__to_test_.xml',
        ...     tag='Bolometer', key='id') # doctest: +NORMALIZE_WHITESPACE
        array(['Cr0Df0Sq0_0t', 'Cr0Df0Sq0_0b', 'Cr0Df0Sq0_1t',
        'Cr0Df0Sq0_1b', 'Cr0Df0Sq0_2t', 'Cr0Df0Sq0_2b',
        'Cr0Df0Sq0_3t', 'Cr0Df0Sq0_3b'], dtype='|S12')

        Return the x coordinates of the 8 bolometers in the focal plane
        >>> fp.unpack_hwmap(fn='./focal_plane__to_test_.xml',
        ...     tag='Bolometer', key='xCoordinate', dtype=float)
        ...     # doctest: +NORMALIZE_WHITESPACE
        array([-15., -15.,  15.,  15., -15., -15.,  15.,  15.])

        """
        tree = ET.parse(fn)
        root = tree.getroot()

        data = []
        for i in range(len(root.getchildren())):
            if root[0].tag == tag:
                data.append(root[i].get(key))
            else:
                for j in range(len(root[0].getchildren())):
                    if root[0][0].tag == tag:
                        data.append(root[i][j].get(key))
                    else:
                        for k in range(len(root[0][0].getchildren())):
                            if root[0][0][0].tag == tag:
                                data.append(root[i][j][k].get(key))
                            else:
                                for l in range(
                                        len(root[0][0][0].getchildren())):
                                    data.append(root[i][j][k][l].get(key))

        return np.array(data, dtype=dtype)

    @staticmethod
    def read_hwmap(fn, tag):
        """
        Routine to show the structure of the hardware map generated
        by create_hwmap(). Particularly useful for interfacing.

        Parameters
        ----------
        fn : string
            xml file containing the hardware map written by create_hwmap.
        tag : string
            Tag corresponding to the layer in the xml file (Crate, DfMuxBoard,
            Squid, Bolometer).

        Returns
        ----------
        data : ndarray
            The data corresponding to `key`.

        Examples
        ----------
        >>> fp = focal_plane(debug=False)

        Return the id of the Crate boards in the focal plane (one here)
        >>> fp.read_hwmap(fn='./focal_plane__to_test_.xml', tag='Crate')
        [['id']]

        Return the id of the DfMux boards in the focal plane (one here)
        >>> fp.read_hwmap(fn='./focal_plane__to_test_.xml', tag='DfMuxBoard')
        ...     # doctest: +NORMALIZE_WHITESPACE
        [['id']]

        Return the id of the Squids in the focal plane (one here)
        >>> fp.read_hwmap(fn='./focal_plane__to_test_.xml', tag='Squid')
        ...     # doctest: +NORMALIZE_WHITESPACE
        [['id']]

        Return the id of the 8 bolometers (4 pairs) in the focal plane
        >>> fp.read_hwmap(fn='./focal_plane__to_test_.xml', tag='Bolometer')
        ...     # doctest: +NORMALIZE_WHITESPACE
        [['xCoordinate', 'focalPlaneIndex', 'yCoordinate',
          'polangle_orientation', 'polarizationMode',
          'polarizationOrientation', 'id', 'channel']]

        """
        tree = ET.parse(fn)
        root = tree.getroot()

        keys = []
        if root[0].tag == tag:
            keys.append(root[0].keys())
        elif root[0][0].tag == tag:
            keys.append(root[0][0].keys())
        elif root[0][0][0].tag == tag:
            keys.append(root[0][0][0].keys())
        elif root[0][0][0][0].tag == tag:
            keys.append(root[0][0][0][0].keys())

        return keys

    @staticmethod
    def show_hwmap(fn_in, fn_out='plot_hardware_map_test.png',
                   save_on_disk=True, display=False):
        """
        Grab the hardware map and show the focal plane of the instrument

        Parameters
        ----------
        fn_in : string
            The harware map to display (xml file).
        fn_out : string, optional
            Name of the output file containing the plot of the focal plane.
            Provide the extension (format: png or pdf).
        save_on_disk : bool
            If True, save the plot on disk.
        display : bool
            If True, show the plot.

        Examples
        ---------
        >>> fp = focal_plane(debug=False)
        >>> fp.show_hwmap(fn_in='./focal_plane__to_test_.xml',
        ...     save_on_disk=False, display=False)
        """
        if not display:
            import matplotlib as mpl
            mpl.use('Agg')
            import matplotlib.pyplot as pl
            pl.ioff()
        else:
            import matplotlib.pyplot as pl

        bolox = focal_plane.unpack_hwmap(
            fn_in, 'Bolometer', 'xCoordinate', dtype=float)
        boloy = focal_plane.unpack_hwmap(
            fn_in, 'Bolometer', 'yCoordinate', dtype=float)
        color = focal_plane.unpack_hwmap(
            fn_in, 'Bolometer', 'polangle_orientation')

        fig, ax = pl.subplots(1, 2, figsize=(10, 5))
        ## Top pixel
        ax[0].scatter(bolox[::2], boloy[::2],
                      c=color[::2], alpha=1, s=30, cmap=pl.cm.jet)
        ax[0].scatter(bolox[::2], boloy[::2],
                      c='black', s=30, marker='|',
                      label='Top pixel', alpha=0.6)
        ax[0].set_ylabel('y position (cm)')
        ax[0].set_xlabel('x position (cm)')
        ax[0].set_title('Top pixels')

        ## Bottom pixel
        ax[1].scatter(bolox[1::2], boloy[1::2],
                      c=color[1::2], alpha=1, s=30, cmap=pl.cm.jet)
        ax[1].scatter(bolox[1::2], boloy[1::2],
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

class beam_model():
    """ Class to handle the beams of the detectors """
    def __init__(self,
                 focal_plane, FWHM=3.5, beam_seed=58347,
                 projected_fp_size=3.,
                 output_folder='./', name='_to_test_', debug=False):
        """
        Parameters
        ----------
        focal_plane : focal_plane instance
            Instance of focal_plane containing focal plane parameters.
        FWHM : float, optional
            Full Width Half Maximum of the beam (in arcmin).
            Default = 3.5 arcmin.
        beam_seed : int
            Seed used to generate angle of rotation of beam axes.
            Default is 58347.
        projected_fp_size : float, optional
            Size of the focal plane on the sky (in degree). This has to
            do with the size of the mirror. Default = 3 degrees.
        output_folder : string, optional
            The folder where the data will be stored.
        name : string, optional
            Tag for the output file (without the extension).
        debug : boolean, optional
            If True, print out a number of useful comments for debugging.
        """
        ## Focal plane parameters
        self.focal_plane = focal_plane

        ## Beam model and mirror parameters
        self.FWHM = FWHM
        self.beam_seed = beam_seed
        self.projected_fp_size = projected_fp_size

        ## Paths and names
        self.output_folder = output_folder
        self.name = name
        self.debug = debug

        ## TODO make output format uniform...
        ## Probably use hdf5 everywhere!
        self.output_file = os.path.join(
            self.output_folder, 'beamprm_' + self.name + '.fits')

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
        >>> fp = focal_plane(debug=False)
        >>> bm = beam_model(fp, debug=False)
        >>> beam = bm.generate_beam_parameters()
        >>> print(beam['xpos']) # doctest: +NORMALIZE_WHITESPACE
        [-0.01308997 -0.01308997 -0.01308997 -0.01308997
          0.01308997  0.01308997  0.01308997  0.01308997
         -0.01308997 -0.01308997 -0.01308997 -0.01308997
          0.01308997  0.01308997  0.01308997  0.01308997]
        """

        beamprm_header = ['Amp', 'Amp_err',
                          'ellip_ang', 'ellip_ang_err',
                          'sig_1', 'sig_1_err',
                          'sig_2', 'sig_2_err',
                          'xpos', 'xpos_err',
                          'ypos', 'ypos_err']
        beamprm_fields = {k: np.zeros(
            self.focal_plane.nbolometer) for k in beamprm_header}

        ## Position of the bolometers
        ## (pairs cm -> bolometers cm -> bolometers radians)
        xp = self.focal_plane.unpack_hwmap(
            self.focal_plane.output_file,
            'Bolometer', 'xCoordinate', dtype=float)
        yp = self.focal_plane.unpack_hwmap(
            self.focal_plane.output_file,
            'Bolometer', 'yCoordinate', dtype=float)

        beamprm_fields['xpos'], beamprm_fields['ypos'] = \
            focal_plane.convert_pair_to_bolometer_position(xp, yp)

        beamprm_fields['xpos'], beamprm_fields['ypos'] = \
            self.convert_cm_to_rad(
            beamprm_fields['xpos'],
            beamprm_fields['ypos'],
            conversion=self.projected_fp_size /
            self.focal_plane.fp_size * np.pi / 180.)

        ## Generate Gaussian beams.
        ## FWHM arcmin -> FWHM rad -> sigma rad
        FWHM_rad = self.FWHM / 60. * np.pi / 180.
        sigma_rad = FWHM_rad / np.sqrt(8 * np.log(2))

        beamprm_fields['sig_1'] = np.ones(
            self.focal_plane.nbolometer) * sigma_rad
        beamprm_fields['sig_2'] = np.ones(
            self.focal_plane.nbolometer) * sigma_rad

        ## Angle of rotation for the ellipses.
        ## Between -90 and 90 degrees.
        state = np.random.RandomState(self.beam_seed)
        beamprm_fields['ellip_ang'] = state.uniform(
            -90, 90, self.focal_plane.nbolometer)

        ## Amplitude of the beams
        ## Default is one.
        beamprm_fields['Amp'] = np.ones(self.focal_plane.nbolometer)

        return beamprm_fields

    @staticmethod
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
        >>> fp = focal_plane(debug=False)
        >>> bm = beam_model(fp, debug=False)
        >>> xcm, ycm = coordinates_on_grid(
        ...                 row_size=fp.fp_size, nx2=fp.npair)
        >>> print(bm.convert_cm_to_rad(xcm, ycm,
        ...     conversion=bm.projected_fp_size / fp.fp_size * np.pi / 180.))
        ...     # doctest: +NORMALIZE_WHITESPACE
        (array([-0.01308997,  0.01308997, -0.01308997,  0.01308997]),
         array([-0.01308997, -0.01308997,  0.01308997,  0.01308997]))

        Note that in this simple example, the focal plane size is 60 cm,
        but the pairs are separated by 30 cm only. For a realistic
        number of pairs (>100), the whole focal plane space is used.
        """
        return np.array(xcm) * conversion, np.array(ycm) * conversion

    @staticmethod
    def construct_beammap(beamprm, ct, cb, nx, pix_size):
        """
        Construct the pixel beam maps
        (sum and difference of bolometer beam maps)

        Parameters
        ----------
        beamprm : dictionary
            Dictionary containing the beam parameters.
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
        >>> fp = focal_plane(debug=False)
        >>> bm = beam_model(fp, debug=False)
        >>> pix_size = 0.5 / 180. * np.pi / 60. # 0.5 arcmin in rad
        >>> summap, diffmap = bm.construct_beammap(
        ...     beamprm=bm.beamprm, ct=0, cb=1, nx=4, pix_size=pix_size)
        >>> print(summap) # doctest: +NORMALIZE_WHITESPACE
        [[ 0.77520676  0.86809111  0.86809111  0.77520676]
         [ 0.86809111  0.97210474  0.97210474  0.86809111]
         [ 0.86809111  0.97210474  0.97210474  0.86809111]
         [ 0.77520676  0.86809111  0.86809111  0.77520676]]
        """
        # Translate beams to origin and maintain differential pointing
        dx = beamprm['xpos'][ct] - beamprm['xpos'][cb]
        dy = beamprm['ypos'][ct] - beamprm['ypos'][cb]

        tx = 0.5 * dx
        bx = -0.5 * dx
        ty = 0.5 * dy
        by = -0.5 * dy

        xy2f = coordinates_on_grid(pix_size=pix_size, nx=nx)

        tmap = beam_model.gauss2d(xy2f, tx, ty,
                                  beamprm['Amp'][ct], beamprm['sig_1'][ct],
                                  beamprm['sig_2'][ct],
                                  beamprm['ellip_ang'][ct]).reshape((nx, nx))

        bmap = beam_model.gauss2d(xy2f, bx, by,
                                  beamprm['Amp'][cb], beamprm['sig_1'][cb],
                                  beamprm['sig_2'][cb],
                                  beamprm['ellip_ang'][cb]).reshape((nx, nx))

        summap = 0.5 * (tmap + bmap)
        diffmap = 0.5 * (tmap - bmap)

        return summap, diffmap

    @staticmethod
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
        >>> fp = focal_plane(debug=False)
        >>> bm = beam_model(fp, debug=False)
        >>> pix_size = 0.5 / 180. * np.pi / 60. # 0.5 arcmin in rad
        >>> xy = coordinates_on_grid(pix_size=pix_size, nx=4)
        >>> beam_model.gauss2d(xy, x_0=0, y_0=0, Amp=1.,
        ...     sig_xp=bm.beamprm['sig_1'][0],
        ...     sig_yp=bm.beamprm['sig_2'][0], psi=0)
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

    def savetodisk(self, headerdict={}):
        """
        Save beam model into disk.
        The data are written into a fits file.

        Parameters
        ----------
        headerdict : dictionary, optional
            Dictionary containing header informations.

        Examples
        ----------
        >>> fp = focal_plane(debug=False)
        >>> bm = beam_model(fp, debug=False)
        >>> bm.savetodisk()
        """

        ## Header of the files
        if len(headerdict) == 0:
            headerdict = {}
            headerdict['ndet'] = self.focal_plane.nbolometer
            headerdict['author'] = 'me'
            headerdict['date'] = str(datetime.date.today())

        savefits_todisk(headerdict, self.beamprm, self.output_file)

class pointing_model():
    """ Class to handle the pointing model of the telescope """
    def __init__(self,
                 pm_name='5params', output_folder='./', name='_to_test_'):
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
        output_folder : string, optional
            The folder where the data will be stored.
        name : string, optional
            Tag for the output file (without the extension).

        Examples
        --------
        >>> pm = pointing_model()
        >>> [(i, round(j,3)) for i, j in zip(
        ...     pm.allowed_params.split(), pm.value_params)]
        ... # doctest: +NORMALIZE_WHITESPACE
        [('ia', -10.285), ('ie', 8.74),
         ('ca', -15.598), ('an', -0.51),
         ('aw', 0.109)]

        >>> pm = pointing_model(pm_name='super-model-trop-bien')
        ... # doctest: +NORMALIZE_WHITESPACE, +ELLIPSIS
        Traceback (most recent call last):
         ...
        ValueError: Only the five-parameter pointing model
        (Mangum 2001) is implemented for the moment (pm_name = 5params)
        """
        self.output_folder = output_folder
        self.name = name
        self.pm_name = pm_name

        ## TODO make output format uniform...
        ## Probably use hdf5 everywhere!
        self.output_file = os.path.join(
            self.output_folder, 'pm_' + self.name + '.fits')

        if self.pm_name == '5params':
            self.five_parameter_pointing_model()
        else:
            raise ValueError('Only the five-parameter ' +
                             'pointing model (Mangum 2001) is implemented ' +
                             'for the moment (pm_name = 5params)')

    def five_parameter_pointing_model(self):
        """
        Parameters based on Polarbear measurements.
        """
        self.allowed_params = 'ia ie ca an aw'

        self.value_params = np.array([-10.28473073, 8.73953334,
                                      -15.59771781, -0.50977716, 0.10858016])

        ## Set this to zero for the moment
        self.RMS_AZ = 0.0
        self.RMS_EL = 0.0
        self.RMS_RESID = 0.0

    def savetodisk(self, headerdict={}):
        """
        Save pointing model into disk.
        The data are written into a fits file.

        Parameters
        ----------
        headerdict : dictionary, optional
            Dictionary containing header informations.

        Examples
        ----------
        >>> pm = pointing_model()
        >>> pm.savetodisk()
        """
        if len(headerdict) == 0:
            headerdict['allowed_params'] = self.allowed_params
            headerdict['RMS AZ'] = str(self.RMS_AZ)
            headerdict['RMS EL'] = str(self.RMS_EL)
            headerdict['RMS RESID'] = str(self.RMS_RESID)

        data = {}
        data['pointing_model'] = self.value_params

        savefits_todisk(headerdict, data, self.output_file)

class polarisation_angle_model():
    """ Class to handle the detector polarisation angle model """
    def __init__(self, focal_plane, output_folder='./', name='_to_test_'):
        """
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
        focal_plane : focal_plane instance
            instance of focal_plane containing focal plane parameters.
        output_folder : string, optional
            The folder where the data will be stored.
        name : string, optional
            Tag for the output file (without the extension).
        """
        self.focal_plane = focal_plane
        self.output_folder = output_folder
        self.name = name
        self.output_file = os.path.join(
            self.output_folder, 'paprm_' + self.name + '.fits')

        self.paprm = self.generate_paprm()

    def generate_paprm(self):
        """
        Construct the polarisation angle parameters: id of the bolometers,
        the polarisation angle for each, and associated errors.
        Errors are set to zero by default.

        Returns
        ----------
        paprm_fields : dictionary
            Dictionary containing the polarisation angle parameters:
            id of the bolometers, the polarisation angle for each,
            and associated errors.

        """
        paprm_header = ['boloid', 'polangle', 'polangle_err']
        paprm_fields = {k: np.zeros(
            self.focal_plane.nbolometer) for k in paprm_header}

        paprm_fields['boloid'] = self.focal_plane.unpack_hwmap(
            fn=self.focal_plane.output_file,
            tag='Bolometer',
            key='id')

        paprm_fields['polangle'] = self.focal_plane.unpack_hwmap(
            fn=self.focal_plane.output_file,
            tag='Bolometer', key='polangle_orientation', dtype=float)

        paprm_fields['polangle_err'] = np.zeros_like(paprm_fields['polangle'])

        return paprm_fields

    def savetodisk(self, headerdict={}):
        """
        Save polarisation angle model into disk.
        The data are written into a fits file.

        Parameters
        ----------
        headerdict : dictionary, optional
            Dictionary containing header informations.

        Examples
        ----------
        >>> fp = focal_plane(debug=False)
        >>> pa = polarisation_angle_model(fp)
        >>> pa.savetodisk()
        """
        ## Header of the files
        if len(headerdict) == 0:
            headerdict = {}
            headerdict['ndet'] = self.focal_plane.nbolometer
            headerdict['author'] = 'me'
            headerdict['date'] = str(datetime.date.today())

        savefits_todisk(headerdict, self.paprm, self.output_file)

class half_wave_plate():
    """ Class to handle the Half-Wave Plate (HWP) """
    def __init__(self, type_HWP='CRHWP', freq_HWP=2., angle_HWP=0.):
        """
        This class provides routines to compute the HWP angles.
        This can be use later to generate timestreams.

        Parameters
        ----------
        type_HWP : string, optional
            The type of HWP that you want to mount on your instrument.
            * CRWHP: continously rotating HWP.
            * stepped: stepped HWP (once a CES).
        freq_HWP : float, optional
            The frequency of rotation of the HWP in Hz.
        angle_HWP : float, optional
            The offset of the HWP in degree.

        """
        self.type_HWP = type_HWP
        self.freq_HWP = freq_HWP
        self.angle_HWP = angle_HWP

        if self.type_HWP not in ['CRHWP', 'stepped']:
            raise ValueError("`type_HWP` has to be 'CRHWP' or 'stepped'.")

        if self.type_HWP is 'stepped' and freq_HWP != 0.0:
            raise AssertionError("You cannot have a stepped HWP and non-" +
                                 "zero frequency! set freq_HWP=0.0 " +
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
        >>> hwp = half_wave_plate(type_HWP='CRHWP', freq_HWP=2., angle_HWP=0.)
        >>> hwp.compute_HWP_angles(sample_rate=100., size=10)
        ... # doctest: +NORMALIZE_WHITESPACE
        array([ 0.        ,  0.12566371,  0.25132741,  0.37699112,  0.50265482,
            0.62831853,  0.75398224,  0.87964594,  1.00530965,  1.13097336])

        Stepped HWP with 30 degrees angle
        >>> hwp = half_wave_plate(type_HWP='stepped',
        ...     freq_HWP=0.0, angle_HWP=30.)
        >>> hwp.compute_HWP_angles(sample_rate=100., size=10)
        ... # doctest: +NORMALIZE_WHITESPACE
        array([ 0.52359878,  0.52359878,  0.52359878,  0.52359878,  0.52359878,
            0.52359878,  0.52359878,  0.52359878,  0.52359878,  0.52359878])

        For a stepped HWP, the frequency must be zero
        >>> hwp = half_wave_plate(type_HWP='stepped',
        ...     freq_HWP=1.0, angle_HWP=30.)
        ... # doctest: +NORMALIZE_WHITESPACE, +ELLIPSIS
        Traceback (most recent call last):
            ...
        AssertionError: You cannot have a stepped HWP and non-zero frequency!
        set freq_HWP=0.0 if you want a stepped HWP.
        """
        angle = self.angle_HWP * np.pi / 180.

        HWP_angles = np.array(
            [angle + t * (self.freq_HWP / sample_rate) *
             2. * np.pi for t in range(size)])

        return HWP_angles

    def update_hardware(self, new_type_HWP, new_freq_HWP, new_angle_HWP):
        """
        Change the behaviour of the HWP.

        Parameters
        ----------
        new_type_HWP : string, optional
            The type of HWP that you want to mount on your instrument.
            * CRWHP: continously rotating HWP.
            * stepped: stepped HWP (once a CES).
        new_freq_HWP : float, optional
            The frequency of rotation of the HWP in Hz.
        new_angle_HWP : float, optional
            The offset of the HWP in degree.

        Examples
        ----------
        Continously rotating HWP at 2 Hz starting at 0 degree.
        >>> hwp = half_wave_plate(type_HWP='CRHWP', freq_HWP=2., angle_HWP=0.)
        >>> hwp.compute_HWP_angles(sample_rate=100., size=10)
        ... # doctest: +NORMALIZE_WHITESPACE
        array([ 0.        ,  0.12566371,  0.25132741,  0.37699112,  0.50265482,
            0.62831853,  0.75398224,  0.87964594,  1.00530965,  1.13097336])

        For some reason, our HWP died!
        >>> hwp.update_hardware(new_type_HWP='stepped',
        ...     new_freq_HWP=0.0, new_angle_HWP=0.0)
        >>> hwp.compute_HWP_angles(sample_rate=100., size=10)
        ... # doctest: +NORMALIZE_WHITESPACE
        array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.])

        """
        self.type_HWP = new_freq_HWP
        self.freq_HWP = new_freq_HWP
        self.angle_HWP = new_angle_HWP


def savefits_todisk(headerdict, datadict, output_file):
    """
    Save data into fits file.

    Parameters
    ----------
    headerdict : dictionary
        Dictionary containing the header of the fits.
    datadict : dictionary
        Dictionary containing the data. Note that the data must be ndarrays.
    output_file : string
        Name of the fits file to create. If it exists already on the disk,
        the previous file is erased before writting the new file.

    Examples
    ----------
    >>> headerdict = {'Take a break': 'read a short story',
    ...               'Author': 'Philip K. Dick'}
    >>> datadict = {'Short Stories': np.array(['Beyond Lies the Wub'])}
    >>> output_file = './erase_me.fits'
    >>> savefits_todisk(headerdict, datadict, output_file)
    >>> os.remove(output_file)
    """
    header_hdu = pyfits.PrimaryHDU()

    # Use HIERARCH.
    for key in headerdict:
        header_hdu.header['HIERARCH {}'.format(key)] = (headerdict[key])

    cols = [pyfits.Column(
                name=label,
                format=dtype_to_fits(datadict[label].dtype),
                array=datadict[label]) for label in datadict]
    data_hdu = pyfits.BinTableHDU.from_columns(pyfits.ColDefs(cols))
    hdulist = pyfits.HDUList([header_hdu, data_hdu])

    if os.path.exists(output_file):
        os.remove(output_file)
    hdulist.writeto(output_file)
    hdulist.close()

def dtype_to_fits(dtype):
    """
    This function will return the appropriate FITS format string from a given
    numpy dtype object.

    Parameters
    ----------
    dtype : type
        Numpy dtype

    Returns
    ----------
    FITS format: str
        The appropriate FITS format string given the numpy dtype object.

    Examples
    ----------
    >>> array = np.ones(2, dtype=int)
    >>> dtype_to_fits(array.dtype)
    'K'
    >>> array = np.ones(2, dtype=str)
    >>> dtype_to_fits(array.dtype)
    'A1'
    """

    # Type dictionary to convert to fits type
    typedict = {'uint8': 'B',
                'uint16': 'I',
                'uint32': 'J',
                'int8': 'A',
                'int16': 'I',
                'int32': 'J',
                'int64': 'K',
                'float32': 'E',
                'float64': 'D',
                'bool': 'L'}

    if (dtype.kind == 'S'):
        # Handle variable-length string types with care:
        return 'A{}'.format(dtype.itemsize)
    else:
        # Use the above list for numeric types.
        return typedict[dtype.name]

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
    remove_test_data(has_id='_to_test_')
