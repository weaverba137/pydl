# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""This module corresponds to the goddard/astro directory in idlutils.
"""
from __future__ import print_function


def airtovac(air):
    """Convert air wavelengths to wavelengths in vacuum.

    Parameters
    ----------
    air : array-like
        Values of wavelength in air in Angstroms.

    Returns
    -------
    array-like
        Values of wavelength in vacuum in Angstroms.

    Notes
    -----
    * Formula from `P. E. Ciddor, Applied Optics, 35, 1566 (1996)
      <http://adsabs.harvard.edu/abs/1996ApOpt..35.1566C>`_.
    * Values of wavelength below 2000 Å are not converted.
    """
    from numpy import zeros
    try:
        vacuum = zeros(air.shape, dtype=air.dtype) + air
        g = vacuum < 2000.0
    except AttributeError:
        # Most likely, vacuum is simply a float.
        vacuum = air
        g = None
        if air < 2000.0:
            return air
    for k in range(2):
        sigma2 = (1.0e4/vacuum)**2
        fact = (1.0 + 5.792105e-2/(238.0185 - sigma2) +
                1.67917e-3/(57.362 - sigma2))
        vacuum = air * fact
    if g is not None:
        vacuum[g] = air[g]
    return vacuum


def gcirc(ra1, dec1, ra2, dec2, units=2):
    """Computes rigorous great circle arc distances.

    Parameters
    ----------
    ra1, dec1, ra2, dec2 : :class:`float` or array-like
        RA and Dec of two points.
    units : { 0, 1, 2 }, optional
        * units = 0: everything is already in radians
        * units = 1: RA in hours, dec in degrees, distance in arcsec.
        * units = 2: RA, dec in degrees, distance in arcsec (default)

    Returns
    -------
    :class:`float` or array-like
        The angular distance.  Units of the value returned depend on the
        input value of `units`.

    Notes
    -----
    The formula below is the one best suited to handling small angular
    separations.  See:
    http://en.wikipedia.org/wiki/Great-circle_distance
    """
    from numpy import arcsin, cos, deg2rad, rad2deg, sin, sqrt
    if units == 0:
        rarad1 = ra1
        dcrad1 = dec1
        rarad2 = ra2
        dcrad2 = dec2
    elif units == 1:
        rarad1 = deg2rad(15.0*ra1)
        dcrad1 = deg2rad(dec1)
        rarad2 = deg2rad(15.0*ra2)
        dcrad2 = deg2rad(dec2)
    elif units == 2:
        rarad1 = deg2rad(ra1)
        dcrad1 = deg2rad(dec1)
        rarad2 = deg2rad(ra2)
        dcrad2 = deg2rad(dec2)
    else:
        raise ValueError('units must be 0, 1 or 2!')
    deldec2 = (dcrad2-dcrad1)/2.0
    delra2 = (rarad2-rarad1)/2.0
    sindis = sqrt(sin(deldec2)*sin(deldec2) +
                  cos(dcrad1)*cos(dcrad2)*sin(delra2)*sin(delra2))
    dis = 2.0*arcsin(sindis)
    if units == 0:
        return dis
    else:
        return rad2deg(dis)*3600.0


def get_juldate(seconds=None):
    """Returns the current Julian date.

    Uses the MJD trick & adds the offset to get JD.

    Parameters
    ----------
    seconds : :class:`int` or :class:`float`, optional
        Time in seconds since the UNIX epoch.  This should only be used
        for testing.

    Returns
    -------
    :class:`float`
        The Julian Day number as a floating point number.

    Notes
    -----
    Do not use this function if high precision is required, or if you are
    concerned about the distinction between UTC & TAI.
    """
    from time import time
    if seconds is None:
        t = time()
    else:
        t = seconds
    mjd = t/86400.0 + 40587.0
    return mjd + 2400000.5


def get_juldate_main():  # pragma: no cover
    """Entry point for the get_juldate command-line script.
    """
    jd = get_juldate()
    print(jd)
    return 0


def vactoair(vacuum):
    """Convert vacuum wavelengths to wavelengths in air.

    Parameters
    ----------
    vacuum : array-like
        Values of wavelength in vacuum in Angstroms.

    Returns
    -------
    array-like
        Values of wavelength in air in Angstroms.

    Notes
    -----
    * Formula from `P. E. Ciddor, Applied Optics, 35, 1566 (1996)
      <http://adsabs.harvard.edu/abs/1996ApOpt..35.1566C>`_.
    * Values of wavelength below 2000 Å are not converted.
    """
    from numpy import zeros
    try:
        air = zeros(vacuum.shape, dtype=vacuum.dtype)
        g = vacuum < 2000.0
    except AttributeError:
        # Most likely, vacuum is simply a float.
        air = vacuum
        g = None
        if vacuum < 2000.0:
            return air
    sigma2 = (1.0e4/vacuum)**2
    fact = (1.0 + 5.792105e-2/(238.0185 - sigma2) +
            1.67917e-3/(57.362 - sigma2))
    air = vacuum/fact
    if g is not None:
        air[g] = vacuum[g]
    return air
