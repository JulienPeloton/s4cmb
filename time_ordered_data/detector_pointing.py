#!/usr/bin/python
"""
Module to simulate pointing of the focal plane detectors.

Reminder mostly for myself:
Az/El specify the direction of a point on the celestial sphere in the
horizontal coordinate system (a spherical coordinate system).
RA/Dec specify the direction of a point on the celestial sphere in the
equatorial coordinate system.

Author: Julien Peloton, j.peloton@sussex.ac.uk
"""
from __future__ import division, absolute_import, print_function

import healpy as hp
import numpy as np
from numpy import cos
from numpy import sin
from numpy import tan
from pyslalib import slalib
import weave

sec2deg = 360.0/86400.0
d2r = np.pi / 180.0
ASTROMETRIC_GEOCENTRIC = 0
APPARENT_GEOCENTRIC = 1
APPARENT_TOPOCENTRIC = 2

class pointing():
    """ Class to handle detector pointing """
    def __init__(self, az_enc, el_enc, time, value_params,
                 allowed_params='ia ie ca an aw',
                 ra_src=0.0, dec_src=0.0, lat=-22.958,
                 ut1utc_fn='data/ut1utc.ephem'):
        """
        Apply pointing model with parameters `value_params` and
        names `allowed_params` to encoder az,el. Order of terms is
        `value_params` is same as order of terms in `allowed_params`.

        Full list of parameters (Thx Fred!):
            an:  azimuth axis tilt north of vertical
            aw:  azimuth axis tilt west of vertical
            an2:  potato chip
            aw2:  potato chip
            npae:  not parallel azimuth/elevation
            ca:  angle of beam to boresight in azimuth
            ia:  azimuth encoder zero
            ie:  elevation encoder zero + angle of beam
                to boresight in elevation
            tf:  cosine flexure
            tfs:  sine flexure
            ref:  refraction
            dt:  timing error in seconds (requires lat argument)
            elt:  time from start elevation correction
            ta1: linear order structural thermal warping in azimuth
            te1: linear order structural thermal warping in elevation
            sa,sa2: solar radiation structural warping in azimuth
            se,se2: solar radiation structural warping in elevation

        Parameters
        ----------
        az_enc : 1d array
            Encoder azimuth in radians.
        el_enc : 1d array
            Encoder elevation in radians.
        time : 1d array
            Encoder time (UTC) in mjd
        value_params : 1d array
            Value of the pointing model parameters (see instrument.py).
            In degrees (see below for full description)
        allowed_params : list of string
            Name of the pointing model parameters used in `value_params`.
        ra_src : float
            RA of the source (center of the patch).
        dec_src : float
            Dec of the source (center of the patch).
        lat : float, optional
            Latitude of the telescope, in degree.
        ut1utc_fn : string
            File containing time correction to UTC.
            The \Delta{UT} (UT1-UTC) is tabulated in IERS circulars
            and elsewhere. It increases by exactly one second at the end
            of each UTC leap second, introduced in order to
            keep \Delta{UT} within \pm 0s.9.
            The 'sidereal \Delta{UT}' which forms part of AOPRMS(13)
            is the same quantity, but converted from solar to sidereal
            seconds and expressed in radians.
            WTF?

        Examples
        ----------
        See hardware.py for more information on the pointing model.
        >>> allowed_params, value_params, az_enc, el_enc, time = \
            pointing.load_fake_pointing()
        >>> pointing = pointing(az_enc, el_enc, time, value_params,
        ...     allowed_params, lat=-22.)
        >>> print(az_enc[2:4], pointing.az[2:4])
        [ 0.12533323  0.18738131] [ 0.11717842  0.17922137]
        """
        self.az_enc = az_enc
        self.el_enc = el_enc
        self.time = time
        self.value_params = value_params
        self.allowed_params = allowed_params
        self.lat = lat * d2r
        self.ut1utc_fn = ut1utc_fn
        self.ra_src = ra_src
        self.dec_src = dec_src

        self.ut1utc = self.get_ut1utc(self.ut1utc_fn, self.time[0])

        ## Initialise the object
        self.az, self.el = self.apply_pointing_model()
        self.azel2radec()

        ## And then for each det, apply offset_detector

    @staticmethod
    def get_ut1utc(ut1utc_fn, mjd):
        """
        Return the time correction to UTC.

        Returns
        ----------
        ut1utc : float
            Contain the time correction to apply to MJD values.

        Examples
        ----------
        >>> round(pointing.get_ut1utc('data/ut1utc.ephem', 56293), 3)
        0.277
        """
        umjds, ut1utcs = np.loadtxt(ut1utc_fn, usecols=(1, 2)).T
        uindex = np.searchsorted(umjds, mjd)
        ut1utc = ut1utcs[uindex]

        return ut1utc

    def apply_pointing_model(self):
        """
        Apply pointing corrections specified by the pointing model.

        Returns
        ----------
        az : 1d array
            The corrected azimuth in arcminutes.
        el : 1d array
            The corrected elevation in arcminutes.
        """
        assert len(self.value_params) == len(self.allowed_params.split()), \
            AssertionError("Vector containing parameters " +
                           "(value_params) has to have the same " +
                           "length than the vector containing names " +
                           "(allowed_params).")

        ## Here are many parameters defining a pointing model.
        ## Of course, we do not use all of them. They are zero by default,
        ## and only those specified by the user will be used.
        params = {p: 0.0 for p in ['an', 'aw', 'an2', 'aw2', 'an4',
                                   'aw4', 'npae', 'ca', 'ia', 'ie', 'tf',
                                   'tfs', 'ref', 'dt', 'elt', 'ta1',
                                   'te1', 'sa', 'se', 'sa2',
                                   'se2', 'sta', 'ste', 'sta2', 'ste2']}

        for param in params:
            if param in self.allowed_params.split():
                index = self.allowed_params.split().index(param)
                params[param] = self.value_params[index]

        params['dt'] *= sec2deg

        ## Azimuth
        azd = -params['an'] * sin(self.az_enc) * sin(self.el_enc)
        azd -= params['aw'] * cos(self.az_enc) * sin(self.el_enc)

        azd -= -params['an2'] * sin(2 * self.az_enc) * sin(self.el_enc)
        azd -= params['aw2'] * cos(2 * self.az_enc) * sin(self.el_enc)

        azd -= -params['an4'] * sin(4 * self.az_enc) * sin(self.el_enc)
        azd -= params['aw4'] * cos(4 * self.az_enc) * sin(self.el_enc)

        azd += params['npae'] * sin(self.el_enc)
        azd -= params['ca']
        azd += params['ia'] * cos(self.el_enc)

        azd += params['dt'] * (
            -sin(self.lat) + cos(self.az_enc) *
            cos(self.lat) * tan(self.el_enc))

        ## Elevation
        eld = params['an'] * cos(self.az_enc)
        eld -= params['aw'] * sin(self.az_enc)
        eld -= params['an2'] * cos(2 * self.az_enc)
        eld -= params['aw2'] * sin(2 * self.az_enc)
        eld -= params['an4'] * cos(4 * self.az_enc)
        eld -= params['aw4'] * sin(4 * self.az_enc)

        eld -= params['ie']
        eld += params['tf'] * cos(self.el_enc)
        eld += params['tfs'] * sin(self.el_enc)
        eld -= params['ref'] / tan(self.el_enc)

        eld += -params['dt'] * cos(self.lat) * sin(self.az_enc)

        eld += params['elt'] * (self.time - np.min(self.time))

        ## Convert back in radian and apply to the encoder values.
        azd *= np.pi / (180.0 * 60.)
        eld *= np.pi / (180.0 * 60.)

        azd /= np.cos(self.el_enc)

        az = self.az_enc - azd
        el = self.el_enc - eld

        return az, el

    def azel2radec(self):
        """
        Given Az/El, time, and time correction returns RA/Dec, parallactic
        angles, and quaternions (used later to rotate detector coordinates).

        Examples
        ----------
        Go from az/el -> ra/dec
        >>> allowed_params, value_params, az_enc, el_enc, time = \
            pointing.load_fake_pointing()
        >>> pointing = pointing(az_enc, el_enc, time, value_params,
        ...     allowed_params, lat=-22.)
        >>> print(pointing.ra[2:4], pointing.dec[2:4])
        [ 3.06982295  3.13757815] [ 0.66201859  0.65269152]
        """
        self.ra, self.dec, self.pa = self.azel2radecpa()
        v_ra = self.ra
        v_dec = self.dec
        v_pa = self.pa
        v_ra_src = self.ra_src
        v_dec_src = self.dec_src

        self.meanpa = np.median(v_pa)

        self.quaternion = quaternion(v_ra, v_dec, v_pa,
                                     v_ra_src, v_dec_src)

        q = self.quaternion.offset_radecpa_makequat()

        assert q.shape == (self.az.size, 4), \
            AssertionError("Wrong size for the quaternions!")

        self.q = q

    def azel2radecpa(self):
        """
        Given Az/El, time, and time correction returns RA/Dec and parallactic
        angles.

        Examples
        ----------
        Go from az/el -> ra/dec/pa
        >>> allowed_params, value_params, az_enc, el_enc, time = \
            pointing.load_fake_pointing()
        >>> pointing = pointing(az_enc, el_enc, time, value_params,
        ...     allowed_params, lat=-22.)
        >>> ra, dec, pa = pointing.azel2radecpa()
        >>> print(ra[2:4], dec[2:4], pa[2:4]) # doctest: +NORMALIZE_WHITESPACE
        [ 3.06982295  3.13757815] [ 0.66201859  0.65269152]
        [-3.00492754 -2.93368246]
        """
        converter = Azel2Radec(self.time[0], self.ut1utc)
        vconv = np.vectorize(converter.azel2radecpa)
        ra, dec, pa = vconv(self.time, self.az, self.el)
        return ra, dec, pa

    def radec2azel(self):
        """
        Given RA/Dec, time, and time correction returns Az/El.

        Examples
        ----------
        Make sure that we go back to our feet az/el -> ra/dec -> az/el
        >>> allowed_params, value_params, az_enc, el_enc, time = \
            pointing.load_fake_pointing()
        >>> pointing = pointing(az_enc, el_enc, time, value_params,
        ...     allowed_params, lat=-22.)
        >>> az, el = pointing.radec2azel()
        >>> assert np.all(np.round(az[2:4],2) == np.round(pointing.az[2:4],2))
        >>> assert np.all(np.round(el[2:4],2) == np.round(pointing.el[2:4],2))
        """
        converter = Azel2Radec(self.time[0], self.ut1utc)
        vconv = np.vectorize(converter.radec2azel)
        az, el = vconv(self.time, self.ra, self.dec)
        return az, el

    def offset_detector(self, azd, eld):
        """
        To compute RA/Dec of each detector from az/el, it is much
        faster to use the quaternions.

        Parameters
        ----------
        azd : 1d array
            The azimuth array for the observation in radian.
        els : 1d array
            The elevation array for the observation in radian.

        Returns
        ----------
        ra : 1d array
            Right ascension in radian.
        dec : 1d array
            Declination in radian.
        pa : 1d array
            Parallactic angle in radian.
        """
        ra, dec, pa = self.quaternion.offset_radecpa_applyquat(
            self.q, -azd, -eld)
        return ra, dec, pa

    @staticmethod
    def load_fake_pointing():
        """
        Load fake pointing parameters for testing purposes.

        Returns
        ----------
        allowed_params : list of string
            Name of the pointing model parameters used in `value_params`.
        value_params : 1d array
            Value of the pointing model parameters (see instrument.py).
            In degrees.
        az_enc : 1d array
            Encoder azimuth in radians.
        el_enc : 1d array
            Encoder elevation in radians.
        time : 1d array
            Encoder time (UTC) in mjd

        """
        allowed_params = 'ia ie ca an aw'
        value_params = [10.28473073, 8.73953334, -15.59771781,
                        -0.50977716, 0.10858016]
        az_enc = np.array([np.sin(2 * np.pi * i / 100)
                           for i in range(100)])
        el_enc = np.ones(100) * 0.5
        time = np.array([56293 + t/84000. for t in range(100)])

        return allowed_params, value_params, az_enc, el_enc, time


class Azel2Radec(object):
    """ Class to handle az/el <-> ra/dec conversion """
    def __init__(self, mjd, ut1utc,
                 lon=-67.786, lat=-22.958, height=5200.,
                 pressure=533.29, temp=273.15, humidity=0.1, epequi=2000.0):
        """
        This class is mostly a wrapper around slalib.
        The default parameters correspond to the Polarbear observation site.

        Parameters
        ----------
        mjd : float
            Date in MJD.
        ut1utc : float
            Time correction (UT1 - UTC) for the date.
        lon : float, optional
            Longitute of the site of observation in degree.
        lat : float, optional
            Latitude of the site of observation in degree.
        height : float, optional
            Altitude of the site of observation.
        pressure : float, optional
            Local atmospheric pressure at the site of observation in millibar.
        temp : float, optional
            Local ambient temperature at the site of observation in K.
        humidity : float, optional
            Local relative humidity at the site of observation
            (value between 0 and 1).
        epequi : float, optional
            Epoch of mean equinox to be used (Julian).

        """
        self.lon = (360. - lon) * d2r
        self.lat = lat * d2r
        self.height = height
        self.pressure = pressure
        self.temp = temp
        self.humidity = humidity
        self.mjd = mjd

        self.epequi = epequi
        self.ut1utc = ut1utc

        self.updateaoprms(mjd)

    def updateaoprms(self, mjd, lapserate=0.0065, xpm=0.0, ypm=0.0):
        """
        Pre-compute the set of apparent to observed place parameters.
        See http://star-www.rl.ac.uk/docs/sun67.htx/sun67ss8.html

        Parameters
        ----------
        lapserate : float
            ta mere
        xpm : float
            polar motion coordinate x (radians)
        ypm : float
            polar motion coordinate y (radians)

        """
        wavelength = (299792458.0 / 150.0e9) * 1e6
        self.aoprms = slalib.sla_aoppa(mjd, self.ut1utc, self.lon, self.lat,
                                       self.height, xpm, ypm, self.temp,
                                       self.pressure, self.humidity,
                                       wavelength, lapserate)

    def azel2radecpa(self, mjd, az, el):
        """
        Given Az/El and time returns RA/Dec and parallactic angle.
        This routine does not return a precisely correct parallactic angle.

        Parameters
        ----------
        mjd : float
            Date in MJD.
        az : float
            Azimuth in radian.
        el : float
            Elevation in radian.

        Returns
        ----------
        ra : float
            Right ascension in radian.
        dec : float
            Declination in radian.
        pa : float
            Parallactic angle in radian.
        """
        zd = np.pi / 2 - el
        amprms = slalib.sla_mappa(self.epequi, mjd)
        self.aoprms = slalib.sla_aoppat(mjd, self.aoprms)

        ra_app1, dec_app1 = slalib.sla_oapqk('a', az, zd + 1e-8, self.aoprms)
        ra1, dec1 = slalib.sla_ampqk(ra_app1, dec_app1, amprms)
        ra_app2, dec_app2 = slalib.sla_oapqk('a', az, zd - 1e-8, self.aoprms)
        ra2, dec2 = slalib.sla_ampqk(ra_app2, dec_app2, amprms)
        pa = slalib.sla_dbear(ra1, dec1, ra2, dec2)
        ra = 0.5 * (ra1 + ra2)
        dec = 0.5 * (dec1 + dec2)

        return ra, dec, pa

    def radec2azel(self, mjd, ra, dec):
        """
        Given RA/Dec and time returns Az/El.

        Parameters
        ----------
        mjd : float
            Date in MJD.
        ra : float
            Right ascension in radian.
        dec : float
            Declination in radian.

        Returns
        ----------
        az : float
            Azimuth in radian.
        el : float
            Elevation in radian.
        """
        amprms = slalib.sla_mappa(self.epequi, mjd)
        self.aoprms = slalib.sla_aoppat(mjd, self.aoprms)
        ra_app, dec_app = slalib.sla_mapqkz(ra, dec, amprms)
        az, zd, a, b, c = slalib.sla_aopqk(ra_app, dec_app, self.aoprms)
        el = np.pi / 2 - zd
        return az, el

class quaternion():
    """ """
    def __init__(self, ra, dec, pa, v_ra_src, v_dec_src):
        """
        """
        self.ra = ra
        self.dec = dec
        self.pa = pa
        self.v_ra_src = v_ra_src
        self.v_dec_src = v_dec_src

    def offset_radecpa_makequat(self):
        """
        """
        qra = self.euler_quatz(self.ra)

        qdec = self.euler_quaty(-self.dec)
        qpa = self.euler_quatx(-self.pa)

        qracen = self.euler_quatz(-self.v_ra_src)
        qdeccen = self.euler_quaty(self.v_dec_src)

        q = self.mult(qdec, qpa)
        q = self.mult(qra, q)
        q = self.mult(qracen, q)
        q = self.mult(qdeccen, q)

        return q

    def offset_radecpa_applyquat(q, azd, eld):
        """
        """
        assert len(q.shape) == 2, AssertionError("Wrong quaternion size!")
        assert q.shape[1] == 4, AssertionError("Wrong quaternion size!")
        qazd = self.euler_quatz(-azd)
        qeld = self.euler_quaty(-eld)

        qpix = self.mult(qazd, qeld)[0]

        # Inlining this is a 30x speed up
        seq = self.mult_inline(q, qpix)

        assert seq.shape[1] == 4, AssertionError("Wrong size!")

        n = seq.shape[0]
        if n > 1024:
            phi, theta, psi = quat_to_radecpa_c(seq)
        else:
            phi, theta, psi = quat_to_radecpa_python(seq)

        return psi, -theta, -phi

    @staticmethod
    def mult(p, q):
        """
        Multiply arrays of quaternions,
        see: http://en.wikipedia.org/wiki/\
            Quaternions#Quaternions_and_the_geometry_of_R3

        Parameters
        ----------
        p : array

        Examples
        ----------
        >>> quaternion.mult(np.array([[3., 4., 5., 2.],
        ...     [2., 2., 2., 2.]]), np.array([1., 2., 3., 4.]))
        ... #doctest +NORMALIZE_WHITESPACE
        array([[ 16.,  16.,  28., -18.],
               [ 12.,   8.,  16.,  -4.]])
        """
        if p.ndim == 1 and q.ndim > 1:
            p = np.tile(p, (q.shape[0], 1))
        if q.ndim == 1 and p.ndim > 1:
            q = np.tile(q, (p.shape[0], 1))
        if q.ndim == 1 and p.ndim == 1:
            p = p.reshape((1, 4))
            q = q.reshape((1, 4))

        ps = p[:, 3]
        qs = q[:, 3]
        pv = p[:, :3]
        qv = q[:, :3]

        pq = np.empty_like(p)
        pq[:, 3] = ps * qs
        pq[:, 3] -= quaternion.arraylist_dot(pv, qv).flatten()

        pq[:, :3] = ps[:, np.newaxis] * qv + \
            pv * qs[:, np.newaxis] + np.cross(pv, qv)

        return pq

    @staticmethod
    def mult_inline(p, q):
        """
        Inline version for when p is an array of quaternions
        and q is a single quaternion. Big speed-up.

        Examples
        ----------
        >>> quaternion.mult_inline(np.array([[3., 4., 5., 2.],
        ...     [2., 2., 2., 2.]]), np.array([1., 2., 3., 4.]))
        ... #doctest +NORMALIZE_WHITESPACE
        array([[ 16.,  16.,  28., -18.],
               [ 12.,   8.,  16.,  -4.]])
        """

        assert p.ndim == 2, AssertionError("Wrong size!")
        assert p.shape[1] == 4, AssertionError("Wrong size!")
        assert q.ndim == 1, AssertionError("Wrong size!")
        assert q.size == 4, AssertionError("Wrong size!")
        pq = np.zeros_like(p)
        n = p.shape[0]

        c_code = '''
        int i;

        for(i=0;i<4*n;i+=4)
        {
            pq[i+3] = p[i+3] * q[3];
            pq[i+3] -= p[i] * q[0] + p[i+1] * q[1] + p[i+2] * q[2];

            pq[i] = p[i+3] * q[0] + p[i] * q[3] + \
                p[i+1] * q[2] - p[i+2] * q[1];

            pq[i+1] = p[i+3] * q[1] + p[i+1] * q[3] + \
                p[i+2] * q[0] - p[i+0] * q[2];

            pq[i+2] = p[i+3] * q[2] + p[i+2] * q[3] + \
                p[i+0] * q[1] - p[i+1] * q[0];
        }
        '''

        shape = pq.shape
        pq = pq.flatten()
        p = p.flatten()
        q = q.flatten()

        weave.inline(c_code, ['pq', 'p', 'q', 'n'])

        pq = pq.reshape(shape)
        return pq

    @staticmethod
    def mult_fortran(p, q):
        """
        Inline version for when p is an array of quaternions
        and q is a single quaternion. Big speed-up.

        Examples
        ----------
        >>> quaternion.mult_fortran(np.array([[3., 4., 5., 2.],
        ...     [2., 2., 2., 2.]]), np.array([1., 2., 3., 4.]))
        ... #doctest +NORMALIZE_WHITESPACE
        array([[ 16.,  16.,  28., -18.],
               [ 12.,   8.,  16.,  -4.]])
        """

        assert p.ndim == 2, AssertionError("Wrong size!")
        assert p.shape[1] == 4, AssertionError("Wrong size!")
        assert q.ndim == 1, AssertionError("Wrong size!")
        assert q.size == 4, AssertionError("Wrong size!")
        pq = np.zeros_like(p)
        n = p.shape[0]

        shape = pq.shape
        pq = pq.flatten()
        p = p.flatten()
        q = q.flatten()

        from detector_pointing_f import detector_pointing_f
        detector_pointing_f.mult_fortran_f(p, q, pq, n)

        pq = pq.reshape(shape)
        return pq

    @staticmethod
    def arraylist_dot(a, b):
        """
        Dot product of ndarrays.
        p and q can have any size.
        Returns a column array.

        Parameters
        ----------
        a : array of arrays
            first list of arrays.
        b : array of arrays
            second list of arrays.

        Examples
        ----------
        >>> a = np.array([[1, 2, 3]])
        >>> b = np.array([[4, 6, 8], [3, 7, 1]])
        >>> quaternion.arraylist_dot(a, b) #doctest: +NORMALIZE_WHITESPACE
        array([[40],
               [20]])
        """
        if a.ndim == 1 and b.ndim == 1:
            return np.array([[np.dot(a, b)]])
        else:
            return np.sum(a * b, axis=1)[:, np.newaxis]

    @staticmethod
    def euler_quatx(alpha):
        """
        Generate quaternion units along x axis

        Parameters
        ----------
        alpha : float
            Polar angle in radian.

        Examples
        ----------
        >>> quaternion.euler_quatx(alpha=np.pi/2.)
        array([ 0.70710678,  0.        ,  0.        ,  0.70710678])

        """
        alpha = np.asarray(alpha)
        z = np.zeros_like(alpha)
        c = np.cos(alpha * 0.5)
        s = np.sin(alpha * 0.5)
        return np.array([s, z, z, c]).T

    @staticmethod
    def euler_quaty(alpha):
        """
        Generate quaternion units along y axis

        Parameters
        ----------
        alpha : float
            Polar angle in radian.

        Examples
        ----------
        >>> quaternion.euler_quaty(alpha=np.pi/2.)
        array([ 0.        ,  0.70710678,  0.        ,  0.70710678])

        """
        alpha = np.asarray(alpha)
        z = np.zeros_like(alpha)
        c = np.cos(alpha * 0.5)
        s = np.sin(alpha * 0.5)
        return np.array([z, s, z, c]).T

    @staticmethod
    def euler_quatz(alpha):
        """
        Generate quaternion units along z axis

        Parameters
        ----------
        alpha : float
            Polar angle in radian.

        Examples
        ----------
        >>> quaternion.euler_quatz(alpha=np.pi/2.)
        array([ 0.        ,  0.        ,  0.70710678,  0.70710678])

        """
        alpha = np.asarray(alpha)
        z = np.zeros_like(alpha)
        c = np.cos(alpha * 0.5)
        s = np.sin(alpha * 0.5)
        return np.array([z, z, s, c]).T


if __name__ == "__main__":
    import doctest
    doctest.testmod()
