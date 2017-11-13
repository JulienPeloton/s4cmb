#!/usr/bin/python
"""
Script to simulate the scan of a CMB experiment.

Author: Julien Peloton, j.peloton@sussex.ac.uk
"""
from __future__ import division, absolute_import, print_function

import os
import ephem
import numpy as np
import healpy as hp
import weave

from s4cmb.scanning_strategy_f import scanning_strategy_f

from pyslalib import slalib

## numerical constants
radToDeg = 180. / np.pi
sidDayToSec = 86164.0905

class ScanningStrategy():
    """ Class to handle the scanning strategy of the telescope """
    def __init__(self, nces=12, start_date='2013/1/1 00:00:00',
                 telescope_longitude='-67:46.816',
                 telescope_latitude='-22:56.396', telescope_elevation=5200.,
                 name_strategy='deep_patch', sampling_freq=30., sky_speed=0.4,
                 ut1utc_fn='s4cmb/data/ut1utc.ephem', language='python',
                 verbose=False):
        """
        A scanning strategy consists in defining the site of observation
        on earth for which we will make the observation, the region
        of the sky to observe, and the schedule of observations.

        Parameters
        ----------
        nces : int, optional
            Number of scans to generate.
        start_date : string, optional
            Starting date for observations. The format is: YYYY/M/D HH:MM:SS.
        telescope_longitude : string, optional
            Longitute (angle) of the telescope. String form: 0:00:00.0.
        telescope_latitude : string, optional
            Latitude (angle) of the telescope. String form: 0:00:00.0.
        telescope_elevation : float, optional
            Height above sea level (in meter).
        name_strategy : string, optional
            Name of a pre-defined scanning strategy to define the boundaries
            of the scan: elevation, azimuth, and time. Only available for the:
            moment name_strategy = deep_patch.
        sampling_freq : float, optional
            Sampling frequency of the bolometers in Hz.
        sky_speed : float, optional
            Azimuth speed of the telescope in deg/s.
        ut1utc_fn : string, optional
            File containing time correction to UTC.
            This is not used here, but pass to the pointing module later on.
        language : string, optional
            Language used for core computations. For big experiments, the
            computational time can be big, and some part of the code can be
            speeded up by interfacing python with C or Fortran.
            Default is python (i.e. no interfacing and can be slow).
            Choose language=C or language=fortran otherwise. Note that C codes
            are compiled on-the-fly (weave), but for fortran codes you need
            first to compile it. See the setup.py or the provided Makefile.
        verbose : bool
            If True, print out several messages to ease the debug.
            Default is False.

        """
        self.nces = nces
        self.start_date = start_date
        self.name_strategy = name_strategy
        self.sampling_freq = sampling_freq
        self.sky_speed = sky_speed
        self.language = language
        self.ut1utc_fn = ut1utc_fn
        self.verbose = verbose

        self.telescope_location = self.define_telescope_location(
            telescope_longitude, telescope_latitude, telescope_elevation)

        self.define_boundary_of_scan()

    def define_telescope_location(self, telescope_longitude='-67:46.816',
                                  telescope_latitude='-22:56.396',
                                  telescope_elevation=5200.):
        """
        Routine to define the site of observation on earth for which
        positions are to be computed. The location of the Polarbear telescope
        is entered as default.

        Parameters
        ----------
        telescope_longitude : str
            Longitute (angle) of the telescope. String form: 0:00:00.0.
        telescope_latitude : str
            Latitude (angle) of the telescope. String form: 0:00:00.0.
        telescope_elevation : float
            Height above sea level (in meter).

        Returns
        ----------
        location : Observer instance
            An `Observer` instance allows you to compute the positions of
            celestial bodies as seen from a particular latitude and longitude
            on the Earth's surface.

        Examples
        ----------
        >>> scan = ScanningStrategy()
        >>> telescope_location = scan.define_telescope_location()
        >>> telescope_location.elevation
        5200.0

        """
        location = ephem.Observer()
        location.long = telescope_longitude
        location.lat = telescope_latitude
        location.elevation = telescope_elevation

        return location

    def define_boundary_of_scan(self):
        """
        Given a pre-defined scanning strategy,
        define the boundaries of the scan: elevation, azimuth, and time.
        For a custom usage (advanced users), modify this routine.

        Examples
        ----------
        >>> scan = ScanningStrategy(name_strategy='deep_patch')
        >>> scan.elevation # doctest: +NORMALIZE_WHITESPACE
        [30.0, 45.5226, 47.7448, 49.967,
         52.1892, 54.4114, 56.6336, 58.8558,
         61.078, 63.3002, 65.5226, 35.2126]

        Only few scanning strategies are currently available:
        >>> scan = ScanningStrategy(name_strategy='toto')
        ... # doctest: +NORMALIZE_WHITESPACE, +ELLIPSIS
        Traceback (most recent call last):
         ...
        ValueError: Only name_strategy = deep_patch or shallow_patch or
        custom are currently available. For another usage (advanced users),
        modify this routine.

        >>> scan.allowed_scanning_strategies
        ['deep_patch', 'shallow_patch', 'custom']
        """
        self.allowed_scanning_strategies = [
            'deep_patch',
            'shallow_patch',
            'custom']

        if self.name_strategy == 'deep_patch':
            self.elevation = [30.0, 45.5226, 47.7448, 49.967,
                              52.1892, 54.4114, 56.6336, 58.8558,
                              61.078, 63.3002, 65.5226, 35.2126]

            self.az_min = [134.2263, 162.3532, 162.3532, 162.3532,
                           162.3532, 162.3532, 162.3532, 162.3532,
                           162.3532, 162.3532, 162.3532, 204.7929]

            self.az_max = [154.2263, 197.3532, 197.3532, 197.3532,
                           197.3532, 197.3532, 197.3532, 197.3532,
                           197.3532, 197.3532, 197.3532, 224.7929]

            self.begin_LST = ['17:07:54.84', '22:00:21.76', '22:00:21.76',
                              '22:00:21.76', '22:00:21.76', '22:00:21.76',
                              '22:00:21.76', '22:00:21.76', '22:00:21.76',
                              '22:00:21.76', '22:00:21.76', '2:01:01.19']

            self.end_LST = ['22:00:21.76', '02:01:01.19', '02:01:01.19',
                            '02:01:01.19', '02:01:01.19', '02:01:01.19',
                            '02:01:01.19', '02:01:01.19', '02:01:01.19',
                            '02:01:01.19', '02:01:01.19', '6:53:29.11']

            self.dec_min = None
            self.dec_max = None
            self.begin_RA = None
            self.end_RA = None
            self.orientation = None

            ## Center of the patch in RA/Dec
            self.ra_mid = 0.
            self.dec_mid = -57.5
        elif self.name_strategy == 'shallow_patch':
            self.elevation = [30., 30., 30., 30., 30., 30.]

            self.dec_min = ["-5:00:00", "-5:00:00", "-15:00:00",
                            "-15:00:00", "-60:00:00", "-60:00:00"]

            self.dec_max = ["20:00:00", "20:00:00", "20:00:00",
                            "20:00:00", "-15:00:00", "-15:00:00"]

            self.begin_RA = ["7:40:00", "7:40:00", "20:00:00",
                             "20:00:00", "20:00:00", "20:00:00"]

            self.end_RA = ["15:20:00", "15:20:00", "5:40:00",
                           "5:40:00", "5:40:00", "5:40:00"]

            self.orientation = ['west', 'east', 'west', 'east', 'west', 'east']

            self.az_min = None
            self.az_max = None
            self.begin_LST = None
            self.end_LST = None

            ## Center of the patch in RA/Dec
            self.ra_mid = 0.
            self.dec_mid = 0.
        elif self.name_strategy == 'custom':
            self.elevation = None
            self.dec_min = None
            self.dec_max = None
            self.begin_RA = None
            self.end_RA = None
            self.orientation = None
            self.az_min = None
            self.az_max = None
            self.begin_LST = None
            self.end_LST = None

            ## Center of the patch in RA/Dec
            self.ra_mid = 0.
            self.dec_mid = 0.
        else:
            raise ValueError("Only name_strategy = deep_patch or " +
                             "shallow_patch or custom are " +
                             "currently available. For another usage " +
                             "(advanced users), modify this routine.")

    def run_one_scan(self, scan_file, scan_number):
        """
        Generate one observation (i.e. one CES) of the telescope.

        Parameters
        ----------
        scan_file : dictionary
            Empty dictionary which will contain the outputs of the scan.
        scan_number : int
            Index of the scan (between 0 and nces - 1).

        Returns
        ----------
        bool : bool
            Returns True if the scan has been generated, and False if the scan
            already exists on the disk.
        """
        ## Check if we have too much/enough information to make a scan
        msg = "You cannot specify azimuth and declination!"
        assert (getattr(self, 'az_min') and not getattr(self, 'dec_min')) or \
            (not getattr(self, 'az_min') and getattr(self, 'dec_min')), msg
        assert (getattr(self, 'az_max') and not getattr(self, 'dec_max')) or \
            (not getattr(self, 'az_max') and getattr(self, 'dec_max')), msg

        msg = "You cannot specify LST and RA!"
        assert (getattr(self, 'begin_LST') and not getattr(self, 'begin_RA')) or \
            (not getattr(self, 'begin_LST') and getattr(self, 'begin_RA')), msg
        assert (getattr(self, 'end_LST') and not getattr(self, 'end_RA')) or \
            (not getattr(self, 'end_LST') and getattr(self, 'end_RA')), msg

        if getattr(self, 'az_min') and getattr(self, 'az_max'):
            msg = 'You need to define timing bounds!'
            assert getattr(self, 'begin_LST') is not None, msg
            assert getattr(self, 'end_LST') is not None, msg
            azLST = True
            RADEC = False
            if self.verbose:
                print("Using azimuth and LST bounds")

        elif getattr(self, 'dec_min') and getattr(self, 'dec_max'):
            msg = 'You need to define RA bounds!'
            assert getattr(self, 'begin_RA') is not None, msg
            assert getattr(self, 'end_RA') is not None, msg
            msg = 'You need to define orientation of scan (east/west)!'
            assert getattr(self, 'orientation') is not None, msg
            azLST = False
            RADEC = True
            if self.verbose:
                print("Using RA and Dec bounds")

        ## Figure out the elevation to run the scan
        el = self.elevation[scan_number]

        ## Define the sampling rate in Hz
        sampling_freq = self.sampling_freq

        #########################################################
        ## Define geometry of the scan
        #########################################################
        if azLST:
            ## Define geometry of the scan by figuring out the azimuth bounds
            az_mean = (
                self.az_min[scan_number] + self.az_max[scan_number]) * 0.5
            az_throw = (self.az_max[scan_number] -
                        self.az_min[scan_number]) / np.cos(el / radToDeg)

        elif RADEC:
            ## If given bounds in declination, make bounds in azimuth
            ## note there is no sanity checking here!
            az_array = np.linspace(0., 180., endpoint=True, num=360.)
            ra_array = np.zeros(az_array.shape)
            dec_array = np.zeros(az_array.shape)
            dec_min = ephem.degrees(self.dec_min[scan_number])
            dec_max = ephem.degrees(self.dec_max[scan_number])
            if(self.orientation[scan_number] == 'west'):
                az_array += 180.

            for i in range(0, az_array.shape[0]):
                ra_array[i], dec_array[i] = \
                    self.telescope_location.radec_of(
                        az_array[i] / radToDeg, el / radToDeg)

            az_allowed = np.asarray(
                [az_array[i] for i in range(0, az_array.shape[0])
                    if (dec_array[i] > dec_min and dec_array[i] < dec_max)])

            if (az_allowed.shape[0] < 2):
                ms = 'Invalid combination of declination bounds and elevation.'
                print(ms)
                exit(1)

            az_max = np.max(az_allowed)
            az_min = np.min(az_allowed)
            az_mean = (az_min + az_max) * 0.5
            az_throw = (az_max - az_min)

        #########################################################
        ## Define the timing bounds!
        #########################################################
        if azLST:
            LST_now = float(
                self.telescope_location.sidereal_time()) / (2 * np.pi)
            begin_LST = float(
                ephem.hours(self.begin_LST[scan_number])) / (2 * np.pi)
            end_LST = float(
                ephem.hours(self.end_LST[scan_number])) / (2 * np.pi)

            if (begin_LST > end_LST):
                begin_LST -= 1.

            ## Reset the date to correspond to the sidereal time to start
            self.telescope_location.date -= (
                (LST_now - begin_LST) * sidDayToSec) * ephem.second

            ## Figure out how long to run the scan for
            num_pts = int((end_LST - begin_LST) * sidDayToSec * sampling_freq)

        if RADEC:
            self.telescope_location.horizon = el * ephem.degree

            # Define a fixed source at the ra, dec that we want to scan
            target_min_ra = ephem.FixedBody()
            target_max_ra = ephem.FixedBody()

            # Figure out where we are looking now
            ra_target, dec_target = self.telescope_location.radec_of(
                az_mean / radToDeg, el / radToDeg)

            # Instantiate the targets
            target_min_ra._dec = dec_target
            target_max_ra._dec = dec_target
            target_min_ra._ra = self.begin_RA[scan_number]
            target_max_ra._ra = self.end_RA[scan_number]

            # Compute initial RA
            target_min_ra.compute(self.telescope_location)
            target_max_ra.compute(self.telescope_location)
            if(self.orientation[scan_number] == 'east'):
                self.telescope_location.date =  \
                    self.telescope_location.next_rising(target_min_ra)

                # Recompute coodinates in the light of change of date
                target_min_ra.compute(self.telescope_location)
                target_max_ra.compute(self.telescope_location)

                ## Update number of time samples for the scan
                num_pts = int(
                    (self.telescope_location.next_rising(
                        target_max_ra) - self.telescope_location.date) /
                    ephem.second * sampling_freq)

            if(self.orientation[scan_number] == 'west'):
                self.telescope_location.date =  \
                    self.telescope_location.next_setting(target_min_ra)

                # Recompute coodinates in the light of change of date
                target_min_ra.compute(self.telescope_location)
                target_max_ra.compute(self.telescope_location)

                ## Update number of time samples for the scan
                num_pts = int(
                    (self.telescope_location.next_setting(
                        target_max_ra) - self.telescope_location.date) /
                    ephem.second * sampling_freq)

        ## Run the scan!
        pb_az_dir = 1.
        upper_az = az_mean + az_throw / 2.
        lower_az = az_mean - az_throw / 2.
        az_speed = self.sky_speed / np.cos(el / radToDeg)
        running_az = az_mean

        ## Initialize arrays
        pb_az_array = np.zeros(num_pts)
        pb_mjd_array = np.zeros(num_pts)
        pb_ra_array = np.zeros(num_pts)
        pb_dec_array = np.zeros(num_pts)
        pb_el_array = np.ones(num_pts) * el

        ## Loop over time samples
        # begin_lst = str(self.telescope_location.sidereal_time())
        # Pad scans 10 seconds on either side
        time_padding = 10.0 / 86400.0

        ## Start of the scan
        pb_az_array[0] = running_az
        pb_mjd_array[0] = date_to_mjd(self.telescope_location.date)

        ## Initialize the time
        scan_file['firstmjd'] = pb_mjd_array[0] - time_padding

        ## Update before starting the loop
        running_az += az_speed * pb_az_dir / sampling_freq
        self.telescope_location.date += ephem.second / sampling_freq

        if self.language == 'python':
            for t in range(1, num_pts):
                ## Set the Azimuth and time
                pb_az_array[t] = running_az

                pb_ra_array[t], pb_dec_array[t] = \
                    self.telescope_location.radec_of(
                    pb_az_array[t] * np.pi / 180.,
                    pb_el_array[t] * np.pi / 180.)

                ## Case to change the direction of the scan
                if(running_az > upper_az):
                    pb_az_dir = -1.
                elif(running_az < lower_az):
                    pb_az_dir = 1.

                running_az += az_speed * pb_az_dir / sampling_freq

                ## Increment the time by one second / sampling rate
                pb_mjd_array[t] = pb_mjd_array[t-1] + \
                    ephem.second / sampling_freq

                ## Increment the time by one second / sampling rate
                self.telescope_location.date += ephem.second / sampling_freq

        elif self.language == 'C':
            c_code = r'''
            int t;
            for (t=1;t<num_pts;t++)
            {
                // Set the Azimuth and time
                pb_az_array[t] = running_az;

                // Case to change the direction of the scan
                if (running_az > upper_az)
                {
                    pb_az_dir = -1.0;
                }
                else if (running_az < lower_az)
                {
                    pb_az_dir = 1.0;
                }

                running_az += az_speed * pb_az_dir / sampling_freq;

                // Increment the time by one second / sampling rate
                pb_mjd_array[t] = pb_mjd_array[t-1] + second / sampling_freq;
            }
            '''
            second = 1./24./3600.
            az_speed = float(az_speed)
            pb_az_dir = float(pb_az_dir)
            sampling_freq = float(sampling_freq)
            running_az = float(running_az)
            upper_az = float(upper_az)
            lower_az = float(lower_az)
            weave.inline(c_code, [
                'num_pts',
                'running_az', 'pb_az_array', 'upper_az',
                'lower_az', 'az_speed', 'pb_az_dir', 'pb_mjd_array',
                'second', 'sampling_freq'], verbose=0)

        elif self.language == 'fortran':
            second = 1./24./3600.
            scanning_strategy_f.run_one_scan_f(
                pb_az_array, pb_mjd_array,
                running_az, upper_az, lower_az, az_speed, pb_az_dir,
                second, sampling_freq, num_pts)

        ## Do not use that for precision - it truncates values
        self.telescope_location.date += num_pts * ephem.second / sampling_freq

        ## Save in file
        scan_file['nces'] = self.nces
        scan_file['CES'] = scan_number
        scan_file['sample_rate'] = sampling_freq
        scan_file['sky_speed'] = self.sky_speed
        scan_file['lastmjd'] = pb_mjd_array[-1] + time_padding

        scan_file['azimuth'] = pb_az_array * np.pi / 180
        scan_file['elevation'] = pb_el_array * np.pi / 180
        scan_file['clock-utc'] = pb_mjd_array

        scan_file['RA'] = pb_ra_array
        scan_file['Dec'] = pb_dec_array

        scan_file['nts'] = len(pb_mjd_array)

        if self.verbose:
            print('+-----------------------------------+')
            print(' CES starts at %s and finishes at %s' % (
                mjd_to_greg(scan_file['firstmjd']),
                mjd_to_greg(scan_file['lastmjd'])))
            print(' It lasts %.3f hours' % (
                (scan_file['lastmjd'] - scan_file['firstmjd']) * 24))
            print('+-----------------------------------+')

        ## Add one day before the next CES (to avoid conflict of time)
        self.telescope_location.date += 24 * ephem.second * 3600

        ## Add the scan into the instance
        # self._update('scan{}'.format(scan_number), scan_file)

        return True

    def run(self):
        """
        Generate all the observations (i.e. all CES) of the telescope.

        Examples
        ----------
        >>> scan = ScanningStrategy(sampling_freq=1., nces=2,
        ...     name_strategy='deep_patch')
        >>> scan.run()
        >>> print(scan.scan0['firstmjd'], scan.scan0['lastmjd'])
        56293.6202546 56293.8230093

        By default, the language used for the core computation is the Python.
        It can be quite slow for heavy configuration, and one can set up
        the language to C or fortran for speeding up the computation (x1000).
        Note that C codes are compiled on-the-fly (weave), but for fortran
        codes you need first to compile it. See the setup.py or
        the provided Makefile.
        >>> scan = ScanningStrategy(sampling_freq=1., nces=2,
        ...     language='fortran', name_strategy='shallow_patch')
        >>> scan.run()
        >>> print(round(scan.scan0['firstmjd'], 2))
        56293.37

        Note that you can create your own scanning strategy. First choose
        the custom ones (set everything to None):
        >>> scan = ScanningStrategy(sampling_freq=1., nces=2,
        ...     language='fortran', name_strategy='custom')

        And then define your own parameters. Example:
        >>> scan.nces = 1
        >>> scan.elevation = [30.]
        >>> scan.az_min = [60.]
        >>> scan.az_max = [100.]
        >>> scan.begin_LST = ['17:00:00.00']
        >>> scan.end_LST = ['22:00:00.00']
        >>> scan.run()

        Note that you To create a scanning strategy within s4cmb,
        you need to specify either
        * Elevations + minimum and maximum azimuths (spatial bounds) +
            begining and end Local Sidereal Times (timing bounds)
        * Elevations + minimum and maximum declinations +
            begining and end Right Ascensions (spatial & timing bounds) +
            orientations (east/west)

        """
        ## Initialise the date and loop over CESes
        self.telescope_location.date = self.start_date
        for CES_position in range(self.nces):
            ## Initialise the starting date of observation
            ## It will be updated then automatically
            setattr(self, 'scan{}'.format(CES_position), {})

            # Create the scan strategy
            self.run_one_scan(
                getattr(self, 'scan{}'.format(CES_position)), CES_position)

    def visualize_my_scan(self, nside, reso=6.9, xsize=900, rot=[0, -57.5],
                          nfid_bolometer=6000, fp_size=180., boost=1.,
                          fullsky=False):
        """
        Simple map-making: project time ordered data into sky maps for
        visualisation. It works only in pure python (i.e. if you set
        language='python' when initialising the scanning_strategy class).

        Parameters
        ----------
        nside : int
            Resolution of the healpix map.
        reso : float
            Resolution of the projected map (pixel size) in arcmin.
        xsize : int
            Number of pixel per row for the projected map.
        rot : 1d array of float
            Coordinates of the center of
            the projected map [lon, lat] in degree
        nfid_bolometer : int
            Fiducial number of bolometers for visualisation.
            By default, the scans are done for one reference bolometer.
        fp_size : float
            Size of the focal plane on the sky, in arcmin.
        boost : int
            boost factor to artificially increase the number of hits.
            It doesn't change the shape of the survey (just the amplitude).

        Outputs
        ----------
            * nhit_loc: 1D array, sky map with cumulative hit counts

        Examples
        ----------
        >>> import matplotlib
        >>> matplotlib.use("Agg") ## remove this if you want to display
        >>> scan = ScanningStrategy(sampling_freq=1., nces=1)
        >>> scan.run()
        >>> scan.visualize_my_scan(512)

        """
        import pylab as pl
        if self.language != 'python':
            raise ValueError("Visualisation is available only in pure " +
                             "python because we do not provide (yet) " +
                             "RA and Dec in C or fortran. Relaunch " +
                             "using language='python' in the class " +
                             "ScanningStrategy.")

        npix = hp.pixelfunc.nside2npix(nside)
        nhit = np.zeros(npix)
        for scan_number in range(self.nces):
            scan = getattr(self, 'scan{}'.format(scan_number))

            num_pts = len(scan['clock-utc'])
            pix_global = hp.pixelfunc.ang2pix(
                nside, (np.pi/2.) - scan['Dec'], scan['RA'])

            ## Boresight pointing healpix maps
            nhit_loc = np.zeros(npix, dtype=np.int32)
            scanning_strategy_f.mapmaking(
                pix_global, nhit_loc, npix, num_pts)

            ## Fake large focal plane with many bolometers for visualisation.
            nhit_loc = convolve_focalplane(nhit_loc, nfid_bolometer,
                                           fp_size, boost)

            nhit += nhit_loc

        if self.verbose:
            print('Stats: nhits = {}/{} (fsky={}%), max hit = {}'.format(
                len(nhit[nhit > 0]),
                len(nhit),
                round(len(nhit[nhit > 0])/len(nhit) * 100, 2),
                int(np.max(nhit))))

        nhit[nhit == 0] = hp.UNSEEN
        if not fullsky:
            hp.gnomview(nhit, rot=rot, reso=reso, xsize=xsize,
                        cmap=pl.cm.viridis,
                        title='nbolos = {}, '.format(nfid_bolometer) +
                        'fp size = {} arcmin, '.format(fp_size) +
                        'nhit boost = {}'.format(boost))
        else:
            hp.mollview(nhit, rot=rot, cmap=pl.cm.viridis,
                        title='nbolos = {}, '.format(nfid_bolometer) +
                        'fp size = {} arcmin, '.format(fp_size) +
                        'nhit boost = {}'.format(boost))
        hp.graticule(verbose=self.verbose)

        pl.show()
        pl.clf()

    def _update(self, name, value):
        """
        Wrapper around setattr function.
        Set a named attribute on an object.

        Parameters
        ----------
        name : string
            The name of the new attribute
        value : *
            The value of the attribute.
        """
        setattr(self, name, value)

def convolve_focalplane(bore_nhits, nbolos, fp_radius_amin, boost):
    """
    Given a hit count map,
    perform the focal plane convolution (that is generate the hit map as
    if there were `nbolos`, and boosted).
    Original author: Neil Goeckner-Wald.
    Modifications by Julien Peloton.

    Parameters
    ----------
    bore_nhits : 1D array
        number of hits per pixel for the reference detector.
    nbolos : int
        total number of bolometers desired.
    fp_radius_amin : float
        radius of the focal plane in arcmin.
    boost : float
        boost factor to artificially increase the number of hits.
        It doesn't change the shape of the survey (just the amplitude).

    Returns
    ----------
    focalplane_nhits : 1D array
        Number of hits for the all the detectors.

    Examples
    ----------
    >>> bore_nhits = np.zeros(hp.nside2npix(128))

    ## Put a patch in the center
    >>> bore_nhits[hp.query_disc(128, hp.ang2vec(np.pi/2, 0.),
    ...     radius=10*np.pi/180.)] = 1.

    ## Increase the number of detector (x100) and make a boost (x10)
    >>> conv_bore_nhits = convolve_focalplane(bore_nhits, nbolos=100,
    ...     fp_radius_amin=180, boost=10)
    >>> print(round(np.max(bore_nhits), 2), round(np.max(conv_bore_nhits), 2))
    1.0 1003.87
    """
    # Now we want to make the focalplane maps
    focalplane_nhits = np.zeros(bore_nhits.shape)

    # Resolution of our healpix map
    nside = hp.npix2nside(focalplane_nhits.shape[0])
    resol_amin = hp.nside2resol(nside, arcmin=True)
    fp_rad_bins = int(fp_radius_amin * 2. / resol_amin)
    fp_diam_bins = (fp_rad_bins * 2) + 1

    # Build the focal plane model and a list of offsets
    (x_fp, y_fp) = np.array(
        np.unravel_index(
            range(0, fp_diam_bins**2),
            (fp_diam_bins, fp_diam_bins))).reshape(
                2, fp_diam_bins, fp_diam_bins) - (fp_rad_bins)
    fp_map = ((x_fp**2 + y_fp**2) < (fp_rad_bins)**2)

    bolo_per_pix = nbolos / float(np.sum(fp_map))

    dRA = np.ndarray.flatten(
        (x_fp[fp_map].astype(float) * fp_radius_amin) / (
            fp_rad_bins * 60. * (180. / (np.pi))))
    dDec = np.ndarray.flatten(
        (y_fp[fp_map].astype(float) * fp_radius_amin) / (
            fp_rad_bins * 60. * (180. / (np.pi))))

    pixels_global = np.array(np.where(bore_nhits != 0)[0], dtype=int)
    for n in pixels_global:
        n = int(n)

        # Compute pointing offsets
        (theta_bore, phi_bore) = hp.pix2ang(nside, n)
        phi = phi_bore + dRA * np.sin(theta_bore)
        theta = theta_bore + dDec

        pixels = hp.ang2pix(nside, theta, phi)
        npix_loc = len(pixels)

        ## Necessary because the values in pixels aren't necessarily unique
        ## This is a poor design choice and should probably be fixed
        scanning_strategy_f.convolve_focalplane_f(
            bore_nhits[n], focalplane_nhits, pixels,
            bolo_per_pix, boost, npix_loc)

    return focalplane_nhits

## Here are a bunch of routines to handle dates...

def date_to_mjd(date):
    """
    Convert date in ephem.date format to MJD.

    Parameters
    ----------
    date : ephem.Date
        Floating point value used by ephem to represent a date.
        The value is the number of days since 1899 December 31 12:00 UT. When
        creating an instance you can pass in a Python datetime instance,
        timetuple, year-month-day triple, or a plain float.
        Run str() on this object to see the UTC date it represents.
        ...
        WTF?

    Returns
    ----------
    mjd : float
        Date in the format MJD.

    Examples
    ----------
    >>> e = ephem.Observer()
    >>> e.date = 0.0 ## 1899 December 31 12:00 UT
    >>> mjd = date_to_mjd(e.date)
    >>> print('DATE={} ->'.format(round(e.date, 2)),
    ...     'MJD={}'.format(round(mjd, 2)))
    DATE=0.0 -> MJD=15019.5
    """
    return greg_to_mjd(date_to_greg(date))

def date_to_greg(date):
    """
    Convert date in ephem.date format to gregorian date.

    Parameters
    ----------
    date : ephem.Date
        Floating point value used by ephem to represent a date.
        The value is the number of days since 1899 December 31 12:00 UT. When
        creating an instance you can pass in a Python datetime instance,
        timetuple, year-month-day triple, or a plain float.
        Run str() on this object to see the UTC date it represents.
        ...
        WTF?

    Returns
    ----------
    greg : string
        Gregorian date (format: YYYYMMDD_HHMMSS)

    Examples
    ----------
    >>> e = ephem.Observer()
    >>> e.date = 0.0 ## 1899 December 31 12:00 UT
    >>> greg = date_to_greg(e.date)
    >>> print('DATE={} ->'.format(round(e.date, 2)),
    ...     'GREG={}'.format(greg))
    DATE=0.0 -> GREG=18991231_120000
    """
    date_ = str(date)
    date_ = str(date.datetime())
    greg = date_.split('.')[0].replace('-',
                                       '').replace(':',
                                                   '').replace(' ',
                                                               '_')

    return greg

def greg_to_mjd(greg):
    """
    Convert gregorian date into MJD.

    Parameters
    ----------
    greg : string
        Gregorian date (format: YYYYMMDD_HHMMSS)

    Returns
    ----------
    mjd : float
        Date in the format MJD.

    Examples
    ----------
    >>> greg = '19881103_000000'
    >>> mjd = greg_to_mjd(greg)
    >>> print('GREG={} ->'.format(greg), 'MJD={}'.format(round(mjd, 2)))
    GREG=19881103_000000 -> MJD=47468.0
    """
    year = int(greg[:4])
    month = int(greg[4:6])
    day = int(greg[6:8])
    hour = int(greg[9:11])
    minute = int(greg[11:13])
    second = int(greg[13:15])

    fracday, status = slalib.sla_dtf2d(hour, minute, second)
    mjd, status = slalib.sla_cldj(year, month, day)
    mjd += fracday

    return mjd

def mjd_to_greg(mjd):
    """
    Convert MJD into gregorian date.

    Parameters
    ----------
    mjd : float
        Date in the format MJD.

    Returns
    ----------
    greg : string
        Gregorian date (format: YYYYMMDD_HHMMSS)

    Examples
    ----------
    >>> mjd = greg_to_mjd('19881103_000000')
    >>> greg = mjd_to_greg(mjd)
    >>> print('MJD={} ->'.format(round(mjd, 2)), 'GREG={}'.format(greg))
    MJD=47468.0 -> GREG=19881103_000000
    """
    year, month, day, fracday, baddate = slalib.sla_djcl(mjd)

    if baddate:
        raise ValueError(BadMJD)

    sign, (hour, minute, second, frac) = slalib.sla_dd2tf(2, fracday)

    s = '{:4d}{:2d}{:2d}_{:2d}{:2d}{:2d}'.format(
        year, month, day, hour, minute, second)
    s = s.replace(' ', '0')

    return s


if __name__ == "__main__":
    import doctest
    doctest.testmod()
