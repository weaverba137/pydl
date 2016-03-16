# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""This module corresponds to the mangle directory in idlutils.
"""
import numpy as np
from collections import namedtuple
from astropy.io import fits
from astropy.extern import six
import astropy.utils as au
from . import PydlutilsException


MangleCaps = namedtuple('MangleCaps', ['X', 'CM'])


class ManglePolygon(object):
    """Simple object to represent a polygon.

    Parameters
    ----------
    row : :class:`~astropy.io.fits.fitsrec.FITS_record`
        A row from a :class:`FITS_polygon` object.

    Attributes
    ----------
    CAPS : :class:`MangleCaps`
        Named tuple containing the ``X`` and ``CM`` attributes.  ``X`` is the
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
        try:
            a0 = args[0]
        except IndexError:
            a0 = None
        if isinstance(a0, fits.fitsrec.FITS_record):
            self.NCAPS = int(args[0]['NCAPS'])
            self.WEIGHT = float(args[0]['WEIGHT'])
            self.PIXEL = int(args[0]['PIXEL'])
            self.STR = float(args[0]['STR'])
            self.USE_CAPS = int(args[0]['USE_CAPS'])
            x = args[0]['XCAPS'][0:self.NCAPS, :].copy()
            assert x.shape == (self.NCAPS, 3)
            cm = args[0]['CMCAPS'][0:self.NCAPS].copy()
            assert cm.shape == (self.NCAPS,)
            self.CAPS = MangleCaps(x, cm)
        elif kwargs:
            if 'x' in kwargs and 'cm' in kwargs:
                xs = kwargs['x'].shape
                cm = kwargs['cm'].shape
                assert xs[0] == cm[0]
            else:
                raise ValueError('Input values are missing!')
            self.NCAPS = xs[0]
            if 'weight' in kwargs:
                self.WEIGHT = float(kwargs['weight'])
            else:
                self.WEIGHT = 1.0
            if 'pixel' in kwargs:
                self.PIXEL = int(kwargs['pixel'])
            else:
                self.PIXEL = -1
            if 'use_caps' in kwargs:
                self.USE_CAPS = int(kwargs['use_caps'])
            else:
                # Use all caps by default
                self.USE_CAPS = (1 << self.NCAPS) - 1
            self.CAPS = MangleCaps(kwargs['x'], kwargs['cm'])
            self.STR = 1.0
        else:
            # Empty polygon object
            self.NCAPS = 0
            self.WEIGHT = 0.0
            self.PIXEL = -1
            self.STR = 0.0
            self.USE_CAPS = 0
            self.CAPS = MangleCaps(None, None)
        return

    @au.lazyproperty
    def cmminf(self):
        """Find the smallest cap in a polygon, accounting for negative
        caps.

        Returns
        -------
        :class:`int`
            The index of the smallest cap.
        """
        cmmin = 2.0
        kmin = -1
        for k in range(self.NCAPS):
            if self.CAPS.CM[k] >= 0:
                cmk = self.CAPS.CM[k]
            else:
                cmk = 2.0 + self.CAPS.CM[k]
            if cmk < cmmin:
                cmmin = cmk
                kmin = k
        return kmin

    def garea(self):
        """Compute the area of a polygon.
        """
        cmmin = self.CAPS.CM[self.cmminf]
        if self.NCAPS >= 2 and cmmin > 1.0:
            np = self.NCAPS + 1
        else:
            np = self.NCAPS
        return cmmin


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
        ``True`` if a cap is used.
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


def read_mangle_polygons(filename):
    """Read a "polygon" format ASCII file in Mangle's own format.  These
    files typically have extension ``.ply`` or ``.pol``.

    Parameters
    ----------
    filename : :class:`str`
        Name of FITS file to read.

    Returns
    -------
    :func:`tuple`
        A tuple containing the list polygons read and a list of strings
        containing the header metadata, if any.
    """
    import re
    with open(filename, universal_newlines=True) as ply:
        data = ply.read()
    lines = data.split("\n")
    try:
        npoly = int(line[0].split()[0])
    except ValueError:
        raise PydlutilsException(("Invalid first line of {0}!  " +
                                  "Are you sure this is a Mangle " +
                                  "polygon file?").format(filename))
    p_lines = [i for i, l in enumerate(lines) if lines.startswith('polygon')]
    header = lines[1:p_lines[0]]
    poly = list()
    r1 = re.compile(r'polygon\s+(\d+)\s+\(([^)]+)\):')
    mtypes = {'str': float, 'weight': float, 'pixel': int, 'caps': int}
    for p in p_lines:
        m = r1.match(lines[p])
        g = m.groups()
        pid = int(g[0])
        meta = g[1].strip().split(',')
        m1 = [m.strip().split()[1] for m in meta]
        m0 = [mtypes[m1[i]](m.strip().split()[0]) for i, m in enumerate(meta)]
        metad = dict(zip(m2,m1))
        metad['x'] = list()
        metad['cm'] = list()
        for cap in lines[p+1:p+1+metad['caps']]:
            data = map(float, re.split(r'\s+', cap.strip()))
            metad['x'].append(map(float, data[0:3]))
            metad['cm'].append(float(data[-1]))
        metad['x'] = np.array(metad['x'])
        assert metad['x'].shape == (metad['caps'], 3)
        metad['cm'] = np.array(metad['cm'])
        assert metad['cm'].shape == (metad['caps'],)
        poly.append(ManglePolygon(**metad))
    return (poly, header)


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
