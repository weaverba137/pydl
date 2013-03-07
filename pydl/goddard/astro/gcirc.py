# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def gcirc(ra1,dec1,ra2,dec2,units=2):
    """Computes rigorous great circle arc distances.

    Parameters
    ----------
    ra1, dec1, ra2, dec2 : float or array_like
        RA and Dec of two points.
    units : { 0, 1, 2 }, optional
        * units = 0: everything is already in radians
        * units = 1: RA in hours, dec in degrees, distance in arcsec.
        * units = 2: RA, dec in degrees, distance in arcsec (default)

    Returns
    -------
    gcirc : float or array_like
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
    delra2 =  (rarad2-rarad1)/2.0
    sindis = sqrt( sin(deldec2)*sin(deldec2) +
        cos(dcrad1)*cos(dcrad2)*sin(delra2)*sin(delra2) )
    dis = 2.0*arcsin(sindis)
    if units == 0:
        return dis
    else:
        return rad2deg(dis)*3600.0

