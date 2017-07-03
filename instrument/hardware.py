#!/usr/bin/python
"""
Script to simulate the focal plane of a CMB experiment.
The hardware configuration is stored in a xml file, and read later on.

Author: Julien Peloton, j.peloton@sussex.ac.uk
"""
from __future__ import division, absolute_import, print_function

import os
import copy
import numpy as np

import xml.etree.ElementTree as ET

class hardware():
    """
    Class to handle the hardware of the instrument:
    focal plane configuration, HWP
    """
    def __init__(self,
                 ncrate, ndfmux_per_crate, nsquid_per_mux, npair_per_squid,
                 fp_size, geometry='square',
                 output_folder='./', tag='test', debug=False):
        """
        Initialise our hardware by creating the focal plane and saving
        the configuration inside a xml file.

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
        fp_size : float
            The size of the focal plane in cm.
        geometry : string, optional
            The shape of the focal plane.
        output_folder : string, optional
            The folder where the xml file will be stored.
        tag : string, optional
            The name of the xml file (without the extension).
        debug : boolean, optional
            If True, print out a number of useful comments for debugging.
        """
        self.ncrate = ncrate
        self.ndfmux_per_crate = ndfmux_per_crate
        self.nsquid_per_mux = nsquid_per_mux
        self.npair_per_squid = npair_per_squid
        self.fp_size = fp_size
        self.geometry = geometry
        self.output_folder = output_folder
        self.tag = tag
        self.output_file = self.tag + '.xml'

        self.debug = debug

        self.create_hwmap()

    def coordinate_in_focal_plane(self, npair):
        """
        Define the position of the pairs in the focal plane
        according to a type of geometry. Being a simple person,
        only squared focal plane are available for the moment.
        For a square focal plane, if the initial number of pairs `npair`
        is not square rootable, then the routine will look for the closest
        square root, and leave holes in place of the extra pairs
        in the focal plane.

        Parameters
        ----------
        npair : int
            The total number of pairs of bolometers in the focal plane.

        Returns
        ----------
        xcoord_reduced : ndarray
            Array of length `npair` containing the coordinate of the pairs
            of bolometers along the x axis
        ycoord_reduced : ndarray
            Array of length `npair` containing the coordinate of the pairs
            of bolometers along the y axis

        Examples
        ----------
        Full focal plane
        >>> h.coordinate_in_focal_plane(npair=4)
        (array([-15., -15.,  15.,  15.]), array([-15.,  15., -15.,  15.]))

        Focal plane with hole
        >>> h.coordinate_in_focal_plane(npair=3)
        (array([-15., -15.,  15.]), array([-15.,  15., -15.]))

        """
        if self.geometry == 'square':
            ## Look for the closest number with square root being an integer
            npair_tmp = copy.copy(npair)
            while True:
                npair_per_row = np.sqrt(npair_tmp)
                if int(npair_per_row) == npair_per_row:
                    npair_per_row = int(npair_per_row)
                    break
                else:
                    npair_tmp += 1

            ## Define x and y coordinates for each pair
            radius_fp = self.fp_size / 2.
            xcoord = [
                -radius_fp + (i * 2 + 1) * radius_fp / npair_per_row for i in
                range(npair_per_row) for j in range(npair_per_row)]
            ycoord = [
                -radius_fp + (i * 2 + 1) * radius_fp / npair_per_row for j in
                range(npair_per_row) for i in range(npair_per_row)]
        else:
            raise ValueError('Wrong geometry! {} '.format(self.geometry) +
                             'is not avalaible. ' +
                             'Only `square` is valid for the moment')

        ## Keep only up the pairs up to the initial number of pairs.
        xcoord_reduced = np.array(xcoord[:npair])
        ycoord_reduced = np.array(ycoord[:npair])

        return xcoord_reduced, ycoord_reduced

    def create_hwmap(self):
        """
        Create the hardware map of the instrument,
        that is the xml file containing the focal plane geometry,
        the bolometers id, the wiring, and so on.
        The terminology used here is taken from the Polarbear experiment.
        The hierarchy is the following:
        * CRATE

            * DFMUX

                * SQUID

                    * BOLOMETER (top & bottom)

        Examples
        ----------
        >>> h.create_hwmap()
        Hardware map generated...
        Hardware map written at ./test.xml
        """
        ## Total number of pairs and bolometers in the focal plane
        npair = self.ncrate * self.ndfmux_per_crate * \
            self.nsquid_per_mux * self.npair_per_squid
        nbolo = npair * 2

        ## Retrieve coordinate of the pairs inside the focal plane
        xcoord, ycoord = self.coordinate_in_focal_plane(npair)

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
                    HWmap[crate][dfmux].set('broadcastPort', '5229')
                    HWmap[crate][dfmux].set('crateSlot', '%d' % crate)
                    HWmap[crate][dfmux].set('isClockMaster', '0')
                    HWmap[crate][dfmux].set('isMulticast', '1')
                    HWmap[crate][dfmux].set(
                        'squidControllerIp', '192.168.0.29')
                    HWmap[crate][dfmux].set('broadcastAddress', '224.0.0.1')
                    HWmap[crate][dfmux].set('ipAddress', '192.168.0.29')
                    HWmap[crate][dfmux].set('id', 'Cr%dDf%d' % (crate, dfmux))
                    HWmap[crate][dfmux].set('revision', '0')

                    for squid in range(self.nsquid_per_mux):
                        ## SQUID
                        ET.SubElement(HWmap[crate][dfmux], 'Squid')
                        HWmap[crate][dfmux][squid].set('flux', '0.0')
                        HWmap[crate][dfmux][squid].set('wire', '1')
                        HWmap[crate][dfmux][squid].set('biasReference', '7.0')
                        HWmap[crate][dfmux][squid].set('biasEnd', '7.5')
                        HWmap[crate][dfmux][squid].set('offset', '3.75')
                        HWmap[crate][dfmux][squid].set('biasStart', '6.0')
                        HWmap[crate][dfmux][squid].set('id', 'Cr%dDf%dSq%d' % (
                            crate, dfmux, squid))

                        for pair in range(self.npair_per_squid):
                            ## BOLOMETER

                            ## Split Q/U
                            boloQ = 2*pair
                            boloU = 2*pair + 1

                            if int(pair_id/np.sqrt(npair)) % 2 == 0:
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
                            HWmap[crate][dfmux][squid][boloQ].set(
                                'lcBoardPad', '64')
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
                            HWmap[crate][dfmux][squid][boloU].set(
                                'lcBoardPad', '64')
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
                                assert bolo_id < nbolo, \
                                    'Hardware map generated...'
                                assert pair_id < npair, \
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
    def unpack_hwmap(fn, tag, key):
        """
        Routine to extract data from the hardware map generated
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

        Returns
        ----------
        data : ndarray
            The data corresponding to `key`.

        Examples
        ----------
        >>> h.create_hwmap()
        Hardware map generated...
        Hardware map written at ./test.xml

        Return the id of the Crate boards in the focal plane (one here)
        >>> h.unpack_hwmap(fn='test.xml', tag='Crate', key='id')
        ['Cr0']

        Return the id of the DfMux boards in the focal plane (one here)
        >>> h.unpack_hwmap(fn='test.xml', tag='DfMuxBoard', key='id')
        ['Cr0Df0']

        Return the id of the Squids in the focal plane (one here)
        >>> h.unpack_hwmap(fn='test.xml', tag='Squid', key='id')
        ['Cr0Df0Sq0']

        Return the id of the 8 bolometers (4 pairs) in the focal plane
        >>> h.unpack_hwmap(fn='test.xml', tag='Bolometer', key='id')
        ...     # doctest: +NORMALIZE_WHITESPACE
        ['Cr0Df0Sq0_0t', 'Cr0Df0Sq0_0b', 'Cr0Df0Sq0_1t', 'Cr0Df0Sq0_1b',
         'Cr0Df0Sq0_2t', 'Cr0Df0Sq0_2b', 'Cr0Df0Sq0_3t', 'Cr0Df0Sq0_3b']

        Return the x coordinates of the 8 bolometers in the focal plane
        >>> h.unpack_hwmap(fn='test.xml', tag='Bolometer', key='xCoordinate')
        ...     # doctest: +NORMALIZE_WHITESPACE
        ['-15.0000', '-15.0000', '-15.0000', '-15.0000',
        '15.0000', '15.0000', '15.0000', '15.0000']

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

        return data

    @staticmethod
    def read_hwmap(fn_in, fn_out='plot_hardware_map_test.png',
                   save_on_disk=True, display=False):
        """
        Grab the Hardware map and plot it

        Parameters
        ----------
        fn_in : string
            The harware map to display (xml file)
        fn_out : string, optional
            Name of the output file containing the plot of the focal plane.
            Provide the extension (format: png or pdf)

        Examples
        ---------
        >>> h.create_hwmap()
        Hardware map generated...
        Hardware map written at ./test.xml
        >>> h.read_hwmap(fn_in='test.xml',
        ...     fn_out='plot_hardware_map_test.png',
        ...     save_on_disk=False, display=False)
        """

        import pylab as pl

        bolox = hardware.unpack_hwmap(fn_in, 'Bolometer', 'xCoordinate')
        boloy = hardware.unpack_hwmap(fn_in, 'Bolometer', 'yCoordinate')
        color = hardware.unpack_hwmap(
            fn_in, 'Bolometer', 'polangle_orientation')

        fig, ax = pl.subplots(1, 2)
        ## Top pixel
        ax[0].scatter(bolox[::2], boloy[::2],
                      c=color[::2], alpha=1, s=30, cmap=pl.cm.jet)
        ax[0].scatter(bolox[::2], boloy[::2],
                      c='black', s=30, marker='|',
                      label='Top pixel', alpha=0.6)
        ax[0].set_ylabel('y position (cm)')
        ax[0].set_xlabel('x position (cm)')
        ax[0].legend()

        ## Bottom pixel
        ax[1].scatter(bolox[1::2], boloy[1::2],
                      c=color[1::2], alpha=1, s=30, cmap=pl.cm.jet)
        ax[1].scatter(bolox[1::2], boloy[1::2],
                      c='black', s=30, marker='_',
                      label='Bottom pixel', alpha=0.6)
        ax[1].set_ylabel('y position (cm)')
        ax[1].set_xlabel('x position (cm)')
        ax[1].legend()

        if save_on_disk:
            pl.savefig(fn_out)
        if display:
            pl.show()
        pl.clf()


if __name__ == "__main__":
    import doctest
    doctest.testmod(
        extraglobs={
            'h': hardware(ncrate=1, ndfmux_per_crate=1,
                          nsquid_per_mux=1, npair_per_squid=4,
                          fp_size=60., geometry='square',
                          output_folder='./', tag='test', debug=True)})
