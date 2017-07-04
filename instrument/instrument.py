#!/usr/bin/python
"""
Script to simulate the hardware of a CMB experiment.
The focal plane configuration is stored in a xml file, and read later on.
The instrument models (pointing, beam, polarisation angle) are stored in
a hdf5 file and read later on.

Author: Julien Peloton, j.peloton@sussex.ac.uk
"""
from __future__ import division, absolute_import, print_function

import os
import copy
import numpy as np

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
    """ Class to simulate the hardware of the instrument """
    def __init__(self,
                 ncrate=1, ndfmux_per_crate=1, nsquid_per_mux=1,
                 npair_per_squid=4, fp_size=60.,
                 FWHM=3.5, beam_seed=58347,
                 projected_fp_size=3.,
                 output_folder='./', name='test', debug=False):
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
        output_folder : string, optional
            The folder where the data will be stored.
        name : string, optional
            Tag for the output file (without the extension).
        debug : boolean, optional
            If True, print out a number of useful comments for debugging.

        """
        self.focal_plane = focal_plane(ncrate, ndfmux_per_crate,
                                       nsquid_per_mux, npair_per_squid,
                                       fp_size, output_folder, name, debug)

        self.beam_parameters = beam_model(self.focal_plane, FWHM, beam_seed,
                                          projected_fp_size, output_folder,
                                          name, debug)

class focal_plane():
    """ Class to handle the focal plane of the instrument. """
    def __init__(self,
                 ncrate=1, ndfmux_per_crate=1, nsquid_per_mux=1,
                 npair_per_squid=4, fp_size=60.,
                 output_folder='./', name='test', debug=False):
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
        self.output_file = 'focal_plane_' + self.name + '.xml'

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
        Hardware map written at ./focal_plane_test.xml
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
                            except AssertionError, e:
                                if self.debug:
                                    print(str(e))
                                max_hit = True
                                break

        tree = ET.ElementTree(HWmap)
        fn = os.path.join(self.output_folder, self.output_file)
        tree.write(fn)
        if self.debug:
            print('Hardware map written at {}'.format(
                os.path.join(self.output_folder, self.output_file)))

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
        >>> fp.unpack_hwmap(fn='./focal_plane_test.xml', tag='Crate', key='id')
        ...     # doctest: +NORMALIZE_WHITESPACE
        array(['Cr0'], dtype='|S3')

        Return the id of the DfMux boards in the focal plane (one here)
        >>> fp.unpack_hwmap(fn='./focal_plane_test.xml',
        ...     tag='DfMuxBoard', key='id') # doctest: +NORMALIZE_WHITESPACE
        array(['Cr0Df0'], dtype='|S6')

        Return the id of the Squids in the focal plane (one here)
        >>> fp.unpack_hwmap(fn='./focal_plane_test.xml',
        ...     tag='Squid', key='id') # doctest: +NORMALIZE_WHITESPACE
        array(['Cr0Df0Sq0'], dtype='|S9')

        Return the id of the 8 bolometers (4 pairs) in the focal plane
        >>> fp.unpack_hwmap(fn='./focal_plane_test.xml',
        ...     tag='Bolometer', key='id') # doctest: +NORMALIZE_WHITESPACE
        array(['Cr0Df0Sq0_0t', 'Cr0Df0Sq0_0b', 'Cr0Df0Sq0_1t',
        'Cr0Df0Sq0_1b', 'Cr0Df0Sq0_2t', 'Cr0Df0Sq0_2b',
        'Cr0Df0Sq0_3t', 'Cr0Df0Sq0_3b'], dtype='|S12')

        Return the x coordinates of the 8 bolometers in the focal plane
        >>> fp.unpack_hwmap(fn='./focal_plane_test.xml',
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
        >>> fp.read_hwmap(fn='./focal_plane_test.xml', tag='Crate')
        [['id']]

        Return the id of the DfMux boards in the focal plane (one here)
        >>> fp.read_hwmap(fn='./focal_plane_test.xml', tag='DfMuxBoard')
        ...     # doctest: +NORMALIZE_WHITESPACE
        [['id']]

        Return the id of the Squids in the focal plane (one here)
        >>> fp.read_hwmap(fn='./focal_plane_test.xml', tag='Squid')
        ...     # doctest: +NORMALIZE_WHITESPACE
        [['id']]

        Return the id of the 8 bolometers (4 pairs) in the focal plane
        >>> fp.read_hwmap(fn='./focal_plane_test.xml', tag='Bolometer')
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
        >>> fp.show_hwmap(fn_in='./focal_plane_test.xml',
        ...     fn_out='plot_hardware_map_test.png',
        ...     save_on_disk=False, display=False)
        """

        import pylab as pl

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
                 output_folder='./', name='test', debug=False):
        """
        Parameters
        ----------
        focal_plane : obj
            Class containing focal plane parameters.
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

        self.beamprm = self.generate_beam_parameters()

    def generate_beam_parameters(self):
        """
        Construct the beam parameters such as position of the centroids,
        ellipticity, angle of rotation, and associated errors.

        Returns
        ----------
        beamprm_fields : dictionnary
            Dictionnary containing beam parameters

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
        beamprm : dictionnary
            Dictionnary containing the beam parameters.
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
        '''
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
        '''

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


if __name__ == "__main__":
    import doctest
    doctest.testmod()
