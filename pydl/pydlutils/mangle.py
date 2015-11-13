# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""This module corresponds to the mangle directory in idlutils.
"""
import numpy as np
from astropy.io import fits
from astropy.extern import six


if six.PY3:
    # This works, but wouldn't `int` intead of `long` be OK to use here?
    # http://python3porting.com/differences.html#long
    long = int


def is_cap_used(use_caps, i):
    """Returns ``True`` if a cap is used.

    Parameters
    ----------
    use_caps : :class:`int`
        Bit mask indicating which cap is used.
    i : :class:`int`
        Number indicating which cap we are interested in.

    Returns
    -------
    is_cap_used : :class:`bool`
    """
    return (use_caps & 1 << i) != 0


def read_fits_polygons(filename):
    """Read a "polygon" format FITS file.

    Unlike the IDL version of this routine, the ``XCAPS`` and ``CMPCAPS``
    columns are *not* replaced with ``poly['CAPS']['X']`` and
    ``poly['CAPS']['CM']``.

    Parameters
    ----------
    filename : :class:`str`
        Name of FITS file to read.

    Returns
    -------
    read_fits_polygons : :class:`~astropy.io.fits.FITS_rec`
        The data contained in HDU 1 of the FITS file.
    """
    with fits.open(filename) as hdulist:
        poly = hdulist[1].data
    # poly['CAPS'] = {'CM': list(), 'X': list()}
    # for k in range(len(poly['NCAPS'])):
    #     poly['CAPS']['CM'].append(poly['CMCAPS'][k,0:poly['NCAPS'][k]])
    #     poly['CAPS']['X'].append(
    #                     poly['XCAPS'][k].reshape((8, 3))[0:poly['NCAPS'][k]])
    # del poly['CMCAPS']
    # del poly['XCAPS']
    return poly


def set_use_caps(x, cm, polygon_use_caps, add=False, tol=1.0e-10,
                 allow_doubles=False, allow_neg_doubles=False):
    """Set the bits in use_caps for a polygon.

    Parameters
    ----------
    x : array-like
        Polygon `x` value.
    cm : array-like
        Polygon `cm` value.
    polygon_use_caps : :class:`int`
        Input value of use_caps.
    add : :class:`bool`, optional
        If ``True``, don't initialize the use_caps value to zero, use the
        value of `polygon_use_caps` instead.
    tol : :class:`float`, optional
        Tolerance used to determine whether two caps are identical.
    allow_doubles : :class:`bool`, optional
        Normally, this routine automatically sets use_caps such that no
        two caps with use_caps set are identical.
    allow_neg_doubles : :class:`bool`, optional
        Normally, two caps that are identical except for the sign of `cm`
        would be set unused.  This inhibits that behaviour.

    Returns
    -------
    set_use_caps : :class:`long`
        Value of use_caps.
    """
    if add:
        use_caps = long(polygon_use_caps)
    else:
        use_caps = long(0)
    t2 = tol**2
    use_caps |= long(2)**len(cm) - long(1)
    if not allow_doubles:
        #
        # Check for doubles
        #
        for i in range(len(cm)):
            if is_cap_used(use_caps, i):
                for j in range(i+1, len(cm)):
                    if is_cap_used(use_caps, j):
                        if np.sum((x[i]-x[j])**2) < t2:
                            if ((np.absolute(cm[i]-cm[j]) < tol) or
                                    ((cm[i] + cm[j]) < tol and not
                                    allow_neg_doubles)):
                                #
                                # Don't use
                                #
                                use_caps -= long(2)**j
    return use_caps
