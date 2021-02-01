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
Script to predict the trajectory of a body (Sun, Moon, ...) on the sky
[(theta, phi) coordinates] over the course of a year.

Author: Julien Peloton, peloton@lal.in2p3.fr
"""
from __future__ import division, absolute_import, print_function

import ephem
from datetime import datetime, date, time, timedelta
import numpy as np


class celestial_trajectory:
    """Predict the trajectory of a body (Sun, Moon, ...) """

    def __init__(
        self,
        body,
        lon_observer,
        lat_observer,
        elevation_observer,
        year,
        tz_offset=0.0,
    ):
        """
        Parameters
        ----------
        body : ephem.<body> instance
            <body> can be Sun(), Moon(), and so on. See ephem.
        lon_observer : str
            Longitute (angle) of the telescope. String form: 0:00:00.0.
        lat_observer : str
            Latitude (angle) of the telescope. String form: 0:00:00.0.
        elevation_observer : float
            Height above sea level (in meter).
        year : int
            Year of observation
        tz_offset : float
            Time correction to UTC. Default is 0.
        """
        self.body = body
        self.lon_observer = lon_observer
        self.lat_observer = lat_observer
        self.elevation_observer = elevation_observer
        self.year = year
        self.tz_offset = tz_offset

        self.months = [
            "Jan",
            "Feb",
            "Mar",
            "Apr",
            "May",
            "Jun",
            "Jul",
            "Aug",
            "Sep",
            "Oct",
            "Nov",
            "Dec",
        ]

        self.altaz = []
        self.radec = []
        self.thetaphi = []

        self.initialise_observer()
        self.traj_of_body_oneyear()

    def initialise_observer(self):
        """
        Set the Longitute and Latitude of the Observer.
        """
        # Define the observer
        self.ob = ephem.Observer()
        self.ob.lat, self.ob.lon = self.lat_observer, self.lon_observer
        self.ob.elevation = self.elevation_observer

    @staticmethod
    def radec2thetaphi(radec):
        """
        Correspondance between RA/Dec and theta/phi coordinate systems.

        Parameters
        ----------
        radec : tuple
            RA and Dec angles in radian.

        Returns
        ----------
        theta : float or 1d array
            Theta angle in radian
        phi : float or 1d array
            Phi angle in radian
        """
        theta = np.pi / 2 - radec[1]
        phi = float(radec[0])
        return (theta, phi)

    def radec_of(self, altaz):
        """
        Get RA/Dec from alt/az.

        Parameters
        ----------
        altaz : tuple
            Tuple containing (alt, az) in radians.

        Returns
        ----------
        RA : float
            RA of the body as seen by the observer [radian].
        Dec : float
            Dec of the body as seen by the observer [radian].
        """
        return self.ob.radec_of(altaz[1], altaz[0])

    def alt_az(self, date):
        """
        Compute altitude (elevation) and azimuth of the body as seen
        by the observer.

        Parameters
        ----------
        date : Python datetime instance
            Floating point value used by ephem to represent a date.
            The value is the number of days since 1899 December 31 12:00 UT.
            When creating an instance you can pass in a Python datetime
            instance, timetuple, year-month-day triple, or a plain float.
            Run str() on this object to see the UTC date it represents.
            ...
            WTF?

        Returns
        ----------
        alt : float
            Altitude (elevation) of the body at the time of observation as
            seen by the observer [radian].
        az : float
            Azimuth of the body at the time of observation as seen
            by the observer [radian].
        """
        self.ob.date = date
        self.body.compute(self.ob)
        return (self.body.alt, self.body.az)

    def traj_of_body_oneday(self, date):
        """
        Compute the coordinate (theta, phi) of the Sun at a particular date

        Parameters
        ----------
        date : Python datetime instance
            Floating point value used by ephem to represent a date.
            The value is the number of days since 1899 December 31 12:00 UT. When
            creating an instance you can pass in a Python datetime instance,
            timetuple, year-month-day triple, or a plain float.
            Run str() on this object to see the UTC date it represents.
            ...
            WTF?

        """
        # Update the date
        date = datetime.combine(date, time(12)) - timedelta(hours=self.tz_offset)
        self.ob.date = date

        # Alt/Az
        ALTAZ = self.alt_az(date)
        self.altaz.append(ALTAZ)

        # RA/Dec
        RADEC = self.radec_of(ALTAZ)
        self.radec.append(RADEC)

        # Theta/Phi
        THETAPHI = self.radec2thetaphi(RADEC)
        self.thetaphi.append(THETAPHI)

    def traj_of_body_oneyear(self):
        """
        Plot the course of the sun over one year.
        We assume each month is made of 28 days.

        Returns
        ----------
        coords : 2D array of floats
            The (theta, phi) coordinate for all days of the year. We compute
            the coordinate for only 28 days per month.
        months : list of strings
            Name of the 12 months.

        Examples
        ----------
        >>> sun = ephem.Sun()
        >>> sun_traj = celestial_trajectory(sun, '-67:46.816',
        ...     '-22:56.396', 5200, 2013)
        >>> print('RA =', round(sun_traj.thetaphi[0][0], 2))
        RA = 1.97

        """
        for month in range(1, 13):
            for day in range(1, 28):
                now = date(self.year, 1, 1).replace(month=month, day=day)
                self.traj_of_body_oneday(now)


if __name__ == "__main__":
    import doctest

    doctest.testmod()
