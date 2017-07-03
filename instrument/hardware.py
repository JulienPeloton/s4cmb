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
from xml.etree.ElementTree import ElementTree

def coordinate_in_focal_plane(npair, fp_size, geometry='square'):
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
    fp_size : float
        The size of the focal plane in cm.
    geometry : string, optional
        The shape of the focal plane.

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
    >>> coordinate_in_focal_plane(npair=4, fp_size=60., geometry='square')
    (array([-15., -15.,  15.,  15.]), array([-15.,  15., -15.,  15.]))

    Focal plane with hole
    >>> coordinate_in_focal_plane(npair=3, fp_size=60., geometry='square')
    (array([-15., -15.,  15.]), array([-15.,  15., -15.]))

    """
    if geometry == 'square':
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
        radius_fp = fp_size / 2.
        xcoord = [
            -radius_fp + (i * 2 + 1) * radius_fp / npair_per_row for i in
            range(npair_per_row) for j in range(npair_per_row)]
        ycoord = [
            -radius_fp + (i * 2 + 1) * radius_fp / npair_per_row for j in
            range(npair_per_row) for i in range(npair_per_row)]
    else:
        raise ValueError('Wrong geometry! {} is not '.format(geometry) +
                         'avalaible. Only `square` is valid for the moment')

    ## Keep only up the pairs up to the initial number of pairs.
    xcoord_reduced = np.array(xcoord[:npair])
    ycoord_reduced = np.array(ycoord[:npair])

    return xcoord_reduced, ycoord_reduced

def create_hwmap(ncrate, ndfmux_per_crate, nsquid_per_mux, npair_per_squid,
                 fp_size, geometry='square',
                 output_folder='./', tag='test', debug=False):
    """
    Create the hardware map of the instrument,
    that is the xml file containing the focal plane geometry,
    the bolometers id, the wiring, and so on.
    The terminology used here is taken from the Polarbear experiment.
    The hierarchy is the following:
    * CRATE

        * DFMUX

            * SQUID

                * BOLOMETER

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

    Examples
    ----------
    Create a focal plane with 1,280 pairs of bolometers.
    >>> create_hwmap(ncrate=4, ndfmux_per_crate=1,
    ...    nsquid_per_mux=1, npair_per_squid=4,
    ...    fp_size=60., geometry='square',
    ...    output_folder='./', tag='test', debug=True)
    Hardware map generated...
    Hardware map written at ./test.xml
    """
    ## Total number of pairs and bolometers in the focal plane
    npair = ncrate * ndfmux_per_crate * nsquid_per_mux * npair_per_squid
    nbolo = npair * 2

    ## Retrieve coordinate of the pairs inside the focal plane
    xcoord, ycoord = coordinate_in_focal_plane(
        npair, fp_size, geometry=geometry)

    HWmap = ET.Element('HardwareMap')

    ## Construct the hardware map
    max_hit = False
    while max_hit is False:
        bolo_id = 0
        pair_id = 0
        polangle = []
        name_bolo = []
        for crate in range(ncrate):
            ## CRATE
            ET.SubElement(HWmap, 'Crate')
            HWmap[crate].set('id', str(crate))

            for dfmux in range(ndfmux_per_crate):
                ## DFMUX
                ET.SubElement(HWmap[crate], 'DfMuxBoard')
                HWmap[crate][dfmux].set('broadcastPort', '5229')
                HWmap[crate][dfmux].set('crateSlot', '%d' % crate)
                HWmap[crate][dfmux].set('isClockMaster', '0')
                HWmap[crate][dfmux].set('isMulticast', '1')
                HWmap[crate][dfmux].set('squidControllerIp', '192.168.0.29')
                HWmap[crate][dfmux].set('broadcastAddress', '224.0.0.1')
                HWmap[crate][dfmux].set('ipAddress', '192.168.0.29')
                HWmap[crate][dfmux].set('id', '%d' % (crate + dfmux))
                HWmap[crate][dfmux].set('revision', '0')

                for squid in range(nsquid_per_mux):
                    ## SQUID
                    ET.SubElement(HWmap[crate][dfmux], 'Squid')
                    HWmap[crate][dfmux][squid].set('flux', '0.0')
                    HWmap[crate][dfmux][squid].set('wire', '1')
                    HWmap[crate][dfmux][squid].set('biasReference', '7.0')
                    HWmap[crate][dfmux][squid].set('biasEnd', '7.5')
                    HWmap[crate][dfmux][squid].set('offset', '3.75')
                    HWmap[crate][dfmux][squid].set('biasStart', '6.0')
                    HWmap[crate][dfmux][squid].set('id', 'C%dD%dSq%d' % (
                        crate, dfmux, squid))

                    for pair in range(npair_per_squid):
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
                        ET.SubElement(HWmap[crate][dfmux][squid], 'Bolometer')
                        HWmap[crate][dfmux][squid][boloQ].set(
                            'focalPlaneIndex', '1')
                        HWmap[crate][dfmux][squid][boloQ].set(
                            'polarizationMode', '1')
                        HWmap[crate][dfmux][squid][boloQ].set(
                            'lcBoardPad', '64')
                        ## Position of the bolometer within the SQUID
                        HWmap[crate][dfmux][squid][boloQ].set(
                            'channel', '%d' % boloQ)

                        name_bolo_tmp = 'C%d_%dt' % (crate, pair_id)
                        name_bolo.append(name_bolo_tmp)
                        HWmap[crate][dfmux][squid][boloQ].set(
                            'polarizationOrientation', '0')
                        HWmap[crate][dfmux][squid][boloQ].set(
                            'id', 'C%d_%dt' % (crate, pair_id))
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
                        ET.SubElement(HWmap[crate][dfmux][squid], 'Bolometer')
                        HWmap[crate][dfmux][squid][boloU].set(
                            'focalPlaneIndex', '1')
                        HWmap[crate][dfmux][squid][boloU].set(
                            'polarizationMode', '1')
                        HWmap[crate][dfmux][squid][boloU].set(
                            'lcBoardPad', '64')
                        ## Position of the bolometer within the SQUID
                        HWmap[crate][dfmux][squid][boloU].set(
                            'channel', '%d' % boloU)

                        name_bolo_tmp = 'C%d_%db' % (crate, pair_id)
                        name_bolo.append(name_bolo_tmp)
                        HWmap[crate][dfmux][squid][boloU].set(
                            'polarizationOrientation', '1')
                        HWmap[crate][dfmux][squid][boloU].set(
                            'id', 'C%d_%db' % (crate, pair_id))
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
                            assert bolo_id < nbolo, 'Hardware map generated...'
                            assert pair_id < npair, 'Hardware map generated...'
                        except AssertionError, e:
                            if debug:
                                print(str(e))
                            max_hit = True
                            break

    tree = ET.ElementTree(HWmap)
    fn = os.path.join(output_folder, tag + '.xml')
    tree.write(fn)
    if debug:
        print('Hardware map written at {}'.format(
            os.path.join(output_folder, tag + '.xml')))

def read_hwmap(fn_in, fn_out='hardware_map_test.png'):
    """
    Grab the Hardware map and plot it

    Parameters
    ----------
    fn_in : string
        The harware map to display (xml file)

    Examples
    ---------
    >>> create_hwmap(ncrate=4, ndfmux_per_crate=1,
    ...    nsquid_per_mux=1, npair_per_squid=16,
    ...    fp_size=60., geometry='square',
    ...    output_folder='./', tag='test', debug=False)
    >>> read_hwmap(fn_in='test.xml', fn_out='hardware_map_test.png')
    """

    import pylab as pl
    tree = ET.parse(fn_in)
    root = tree.getroot()

    bolox, boloy, color = [], [], []
    for i in range(len(root.getchildren())):
        for j in range(len(root[0].getchildren())):
            for k in range(len(root[0][0].getchildren())):
                for l in range(len(root[0][0][0].getchildren())):
                    bolox.append(root[i][j][k][l].get('xCoordinate'))
                    boloy.append(root[i][j][k][l].get('yCoordinate'))
                    color.append(root[i][j][k][l].get('polangle_orientation'))

    fig, ax = pl.subplots(1, 2)
    ## Top pixel
    ax[0].scatter(bolox[::2], boloy[::2],
                  c=color[::2], alpha=1, s=30, cmap=pl.cm.jet)
    ax[0].scatter(bolox[::2], boloy[::2],
                  c='black', s=30, marker='|', label='Top pixel', alpha=0.6)
    ax[0].set_ylabel('y position (cm)')
    ax[0].set_xlabel('x position (cm)')
    ax[0].legend()

    ## Bottom pixel
    ax[1].scatter(bolox[1::2], boloy[1::2],
                  c=color[1::2], alpha=1, s=30, cmap=pl.cm.jet)
    ax[1].scatter(bolox[1::2], boloy[1::2],
                  c='black', s=30, label='Bottom pixel', marker='_', alpha=0.6)
    ax[1].set_ylabel('y position (cm)')
    ax[1].set_xlabel('x position (cm)')
    ax[1].legend()
    pl.savefig(fn_out)
    pl.clf()


if __name__ == "__main__":
    import doctest
    doctest.testmod()
