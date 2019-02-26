# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""This module corresponds to the goddard/astro directory in idlutils.
"""
from time import time
import numpy as np
from astropy.units import Angstrom


def airtovac(air):
    """Convert air wavelengths to wavelengths in vacuum.

    Parameters
    ----------
    air : array-like
        Values of wavelength in air in Angstroms.
        :class:`~astropy.units.Quantity` objects with valid length
        dimensions will be internally converted to Angstrom.

    Returns
    -------
    array-like
        Values of wavelength in vacuum in Angstroms.  If a
        :class:`~astropy.units.Quantity` object was passed in, the output
        will be converted to the same units as the input.

    Notes
    -----
    * Formula from `P. E. Ciddor, Applied Optics, 35, 1566 (1996)
      <http://adsabs.harvard.edu/abs/1996ApOpt..35.1566C>`_.
    * Values of wavelength below 2000 Å are not converted.
    """
    try:
        u = air.unit
    except AttributeError:
        u = None
    try:
        t = air.dtype
    except AttributeError:
        # Most likely, air is simply a float.
        t = None
    if t is None:
        if air < 2000.0:
            return air
        vacuum = air
        a = air
        g = None
    else:
        try:
            a = air.to(Angstrom).value
        except AttributeError:
            a = air
        g = a < 2000.0
        if g.all():
            return air
        vacuum = np.zeros(air.shape, dtype=t) + a
    for k in range(2):
        sigma2 = (1.0e4/vacuum)**2
        fact = (1.0 + 5.792105e-2/(238.0185 - sigma2) +
                1.67917e-3/(57.362 - sigma2))
        vacuum = a * fact
    if g is not None:
        vacuum[g] = a[g]
    if u is not None:
        vacuum = (vacuum * Angstrom).to(u)
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
    if units == 0:
        rarad1 = ra1
        dcrad1 = dec1
        rarad2 = ra2
        dcrad2 = dec2
    elif units == 1:
        rarad1 = np.deg2rad(15.0*ra1)
        dcrad1 = np.deg2rad(dec1)
        rarad2 = np.deg2rad(15.0*ra2)
        dcrad2 = np.deg2rad(dec2)
    elif units == 2:
        rarad1 = np.deg2rad(ra1)
        dcrad1 = np.deg2rad(dec1)
        rarad2 = np.deg2rad(ra2)
        dcrad2 = np.deg2rad(dec2)
    else:
        raise ValueError('units must be 0, 1 or 2!')
    deldec2 = (dcrad2-dcrad1)/2.0
    delra2 = (rarad2-rarad1)/2.0
    sindis = np.sqrt(np.sin(deldec2)*np.sin(deldec2) +
                  np.cos(dcrad1)*np.cos(dcrad2)*np.sin(delra2)*np.sin(delra2))
    dis = 2.0*np.arcsin(sindis)
    if units == 0:
        return dis
    else:
        return np.rad2deg(dis)*3600.0


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
    if seconds is None:
        t = time()
    else:
        t = seconds
    mjd = t/86400.0 + 40587.0
    return mjd + 2400000.5


def get_juldate_main():  # pragma: no cover
    """Entry point for the get_juldate command-line script.
    """
    print(get_juldate())
    return 0


def vactoair(vacuum):
    """Convert vacuum wavelengths to wavelengths in air.

    Parameters
    ----------
    vacuum : array-like
        Values of wavelength in vacuum in Angstroms.
        :class:`~astropy.units.Quantity` objects with valid length
        dimensions will be internally converted to Angstrom.

    Returns
    -------
    array-like
        Values of wavelength in air in Angstroms.
        :class:`~astropy.units.Quantity` object was passed in, the output
        will be converted to the same units as the input.

    Notes
    -----
    * Formula from `P. E. Ciddor, Applied Optics, 35, 1566 (1996)
      <http://adsabs.harvard.edu/abs/1996ApOpt..35.1566C>`_.
    * Values of wavelength below 2000 Å are not converted.
    """
    try:
        u = vacuum.unit
    except AttributeError:
        u = None
    try:
        t = vacuum.dtype
    except AttributeError:
        # Most likely, vacuum is simply a float.
        t = None
    if t is None:
        if vacuum < 2000.0:
            return vacuum
        air = vacuum
        v = vacuum
        g = None
    else:
        try:
            v = vacuum.to(Angstrom).value
        except AttributeError:
            v = vacuum
        g = v < 2000.0
        if g.all():
            return vacuum
        air = np.zeros(vacuum.shape, dtype=t) + v
    sigma2 = (1.0e4/v)**2
    fact = (1.0 + 5.792105e-2/(238.0185 - sigma2) +
            1.67917e-3/(57.362 - sigma2))
    air = v / fact
    if g is not None:
        air[g] = v[g]
    if u is not None:
        air = (air * Angstrom).to(u)
    return air
