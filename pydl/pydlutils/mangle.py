# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""This module corresponds to the mangle directory in idlutils.
"""
import numpy as np
from astropy.io import fits
from astropy.extern import six


class ManglePolygon(object):
    """Simple object to represent a polygon.

    Parameters
    ----------
    row : :class:`astropy.iofits.fitsrec.FITS_record`
        A row from a :class:`FITS_polygon` object.

    Attributes
    ----------
    CAPS : :class:`numpy.recarray`
        Array containing the ``X`` and ``CM`` attributes.  ``X`` is the
        direction of the cap on the unit sphere, and ``CM`` is the
        cap's size.
    NCAPS : :class:`int`
        Number of caps in the polygon.
    PIXEL : :class:`int`
        Pixel this polygon is in.
    STR : :class:`float`
        Solid angle of this polygon (steradians).
    USE_CAPS : :class:`int`
        Bitmask indicating which caps to use.
    WEIGHT : :class:`float`
        Weight factor assigned to the polygon.
    """

    def __init__(self, *args, **kwargs):
        if isinstance(args[0], fits.fitsrec.FITS_record):
            self.NCAPS = int(args[0]['NCAPS'])
            self.WEIGHT = float(args[0]['WEIGHT'])
            self.PIXEL = int(args[0]['PIXEL'])
            self.STR = float(args[0]['STR'])
            self.USE_CAPS = int(args[0]['USE_CAPS'])
            x = args[0]['XCAPS'][0:self.NCAPS, :].copy()
            cm = args[0]['CMCAPS'][0:self.NCAPS].copy()
            caps = np.empty((1,), dtype=[('X', x.dtype, x.shape),
                                         ('CM', cm.dtype, cm.shape)])
            self.CAPS = caps.view(np.recarray)
            self.CAPS.X = x
            self.CAPS.CM = cm
        return


class FITS_polygon(fits.FITS_rec):
    """Handle polygons read in from a FITS file.

    This class creates a copy of the ``XCAPS`` and ``CMCAPS`` columns.
    This isn't especially efficient, but we can optimize later.
    """

    #
    # Right now, this class is only instantiated by calling .view() on
    # a FITS_rec object, so only __array_finalize__ is needed.
    #
    # def __new__(*args):
    #     self = super(FITS_polygon, self).__new__(*args)
    #     self._caps = None
    #     return self

    def __array_finalize__(self, obj):
        super(FITS_polygon, self).__array_finalize__(obj)
        self._caps = None

    def __getitem__(self, key):
        if isinstance(key, six.string_types):
            if key.upper() == 'CAPS':
                #
                # 'CAPS' is a container with two columns: 'X', 'CM'.
                #
                if self._caps is None:
                    xdt = self['XCAPS'].dtype
                    cmdt = self['CMCAPS'].dtype
                    self._caps = np.empty((self.size,),
                                          dtype=[('X', xdt, (9, 3)),
                                                 ('CM', cmdt, (9,))]
                                         ).view(np.recarray)
                    self._caps['X'] = self['XCAPS']
                    self._caps['CM'] = self['CMCAPS']
                return self._caps
            return super(FITS_polygon, self).__getitem__(key)
        else:
            return super(FITS_polygon, self).__getitem__(key)

    def __getattribute__(self, key):
        if key.upper() == 'CAPS':
            return self.__getitem__(key)
        return super(FITS_polygon, self).__getattribute__(key)


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
    :class:`bool`
    """
    return (use_caps & 1 << i) != 0


def read_fits_polygons(filename, convert=False):
    """Read a "polygon" format FITS file.

    This function returns a subclass of :class:`~astropy.io.fits.FITS_rec`
    that simulates a "subcolumn" for compatibility with IDL code.
    For example, if ``poly`` is the object returned by this function, then
    the column ``XCAPS`` is accessible as ``poly.CAPS.X``.

    Parameters
    ----------
    filename : :class:`str`
        Name of FITS file to read.
    convert : :class:`bool`, optional
        If ``True``, convert the data to a list of :class:`ManglePolygon`
        objects.

    Returns
    -------
    :class:`~pydl.pydlutils.mangle.FITS_polygon` or :class:`list`
        The data contained in HDU 1 of the FITS file.
    """
    with fits.open(filename, uint=True) as hdulist:
        data = hdulist[1].data
    if convert:
        poly = list()
        for k in range(data.size):
            poly.append(ManglePolygon(data[k]))
    else:
        poly = data.view(FITS_polygon)
    return poly


def set_use_caps(polygon, index_list, add=False, tol=1.0e-10,
                 allow_doubles=False, allow_neg_doubles=False):
    """Set the bits in USE_CAPS for a set of polygons.

    Parameters
    ----------
    polygon : polygon-like
        A polygon object, such as :class:`~pydl.pydlutils.mangle.FITS_polygon`.
    index_list : array-like
        A list of indices of caps to set in the polygon.  Should be no
        longer, nor contain indices greater than the number of caps
        (``NCAPS``).
    add : :class:`bool`, optional
        If ``True``, don't initialize the use_caps value to zero, use the
        existing value associated with `polygon` instead.
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
    :class:`int`
        Value of use_caps.
    """
    if not add:
        polygon.USE_CAPS = 0
    t2 = tol**2
    for i in index_list:
        polygon.USE_CAPS |= (1 << index_list[i])
    if not allow_doubles:
        #
        # Check for doubles
        #
        for i in range(polygon.NCAPS):
            if is_cap_used(polygon.USE_CAPS, i):
                for j in range(i+1, polygon.NCAPS):
                    if is_cap_used(polygon.USE_CAPS, j):
                        if np.sum((polygon.CAPS.X[i, :] -
                                   polygon.CAPS.X[j, :])**2) < t2:
                            if ((np.absolute(polygon.CAPS.CM[i] -
                                             polygon.CAPS.CM[j]) < tol) or
                                ((polygon.CAPS.CM[i] +
                                  polygon.CAPS.CM[j]) < tol and
                                  not allow_neg_doubles)):
                                #
                                # Don't use
                                #
                                polygon.USE_CAPS -= 1 << j
    return polygon.USE_CAPS
