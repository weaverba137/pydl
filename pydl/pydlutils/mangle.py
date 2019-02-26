# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""This module corresponds to the mangle directory in idlutils.

Mangle_ [1]_ is a software suite that supports the concept of polygons on
a sphere, and is used to, for example, describe the `window function`_ of the
Sloan Digital Sky Survey.

This implementation is intended to support the portions of Mangle that
are included in idlutils.  To simplify the coding somewhat, unlike
idlutils, the caps information is accessed through ``polygon.x`` and
``polygon.cm``, not ``polygon.caps.x`` or ``polygon.caps.cm``.

Window function operations are already supported by :func:`~is_in_polygon` and
:func:`~is_in_window`.  However, calculation of polygon area (solid angle)
from first principles is not yet implemented.

Note that in traditional geometry "spherical polygon" means a figure
bounded by *great circles*.  Mangle allows polygons to be bounded by
*any* circle, great or not.  Use care when comparing formulas in this
module to formulas in the mathematical literature.

.. _Mangle: http://space.mit.edu/~molly/mangle/
.. _`window function`: https://www.sdss.org/dr14/algorithms/resolve/

References
----------

..  [1] `Swanson, M. E. C.; Tegmark, Max; Hamilton, Andrew J. S.;
    Hill, J. Colin, 2008 MNRAS 387, 1391
    <http://adsabs.harvard.edu/abs/2008MNRAS.387.1391S>`_.
"""
import re
import warnings
import numpy as np
from astropy.io import fits
# import astropy.utils as au
from . import PydlutilsException, PydlutilsUserWarning


class PolygonList(list):
    """A :class:`list` that contains :class:`ManglePolygon` objects and
    possibly some metadata.

    Parameters
    ----------
    header : :class:`list`, optional
        Set the `header` attribute.

    Attributes
    ----------
    header : :class:`list`
        A list of strings containing metadata.
    """
    def __init__(self, *args, **kwargs):
        super(PolygonList, self).__init__(*args)
        if 'header' in kwargs:
            self.header = kwargs['header']
        else:
            self.header = list()
        return


class ManglePolygon(object):
    """Simple object to represent a polygon.

    A polygon may be instantiated with a row
    (:class:`~astropy.io.fits.fitsrec.FITS_record`) from a
    :class:`FITS_polygon` object, another :class:`ManglePolygon` object
    (copy constructor), keyword arguments, or with no arguments at all,
    in which case it represents the whole sky.

    Attributes
    ----------
    cm : :class:`~numpy.ndarray`
        The size of each cap in the polygon.
    id : :class:`int`
        An arbitrary ID number.
    pixel : :class:`int`
        Pixel this polygon is in.
    use_caps : :class:`int`
        Bitmask indicating which caps to use.
    weight : :class:`float`
        Weight factor assigned to the polygon.
    """

    def __init__(self, *args, **kwargs):
        try:
            a0 = args[0]
        except IndexError:
            a0 = None
        if isinstance(a0, fits.fitsrec.FITS_record):
            self._ncaps = int(a0['NCAPS'])
            self.weight = float(a0['WEIGHT'])
            self.pixel = int(a0['PIXEL'])
            try:
                self.id = int(a0['IFIELD'])
            except KeyError:
                self.id = -1
            self._str = float(a0['STR'])
            self.use_caps = int(a0['USE_CAPS'])
            xcaps = a0['XCAPS']
            if xcaps.shape == (3, ):
                self._x = np.zeros((1, 3), dtype=xcaps.dtype) + xcaps
            else:
                self._x = xcaps[0:self._ncaps, :].copy()
            # assert self._x.shape == (self._ncaps, 3)
            cmcaps = a0['CMCAPS']
            if cmcaps.shape == ():
                self.cm = np.zeros((1,), dtype=cmcaps.dtype) + cmcaps
            else:
                self.cm = cmcaps[0:self._ncaps].copy()
            # assert self.cm.shape == (self._ncaps, )
        elif isinstance(a0, ManglePolygon):
            self._ncaps = a0._ncaps
            self.weight = a0.weight
            self.pixel = a0.pixel
            self.id = a0.id
            self._str = a0._str
            self.use_caps = a0.use_caps
            self._x = a0._x.copy()
            self.cm = a0.cm.copy()
        elif kwargs:
            if 'x' in kwargs and 'cm' in kwargs:
                xs = kwargs['x'].shape
                cm = kwargs['cm'].shape
                assert xs[0] == cm[0]
            else:
                raise ValueError('Input values are missing!')
            self._ncaps = xs[0]
            if 'weight' in kwargs:
                self.weight = float(kwargs['weight'])
            else:
                self.weight = 1.0
            if 'pixel' in kwargs:
                self.pixel = int(kwargs['pixel'])
            else:
                self.pixel = -1
            if 'id' in kwargs:
                self.id = int(kwargs['id'])
            else:
                self.id = -1
            if 'use_caps' in kwargs:
                self.use_caps = int(kwargs['use_caps'])
            else:
                # Use all caps by default
                self.use_caps = (1 << self.ncaps) - 1
            self._x = kwargs['x'].copy()
            self.cm = kwargs['cm'].copy()
            if 'str' in kwargs:
                self._str = float(kwargs['str'])
            else:
                self._str = None
        else:
            #
            # An "empty" polygon represents the whole sky.
            #
            self._ncaps = 0
            self.weight = 1.0
            self.pixel = -1
            self.id = -1
            self.use_caps = 0
            self._x = None
            self.cm = None
            self._str = 4.0*np.pi
        return

    @property
    def ncaps(self):
        """Number of caps in the polygon.
        """
        return self._ncaps

    @property
    def x(self):
        """The orientation of each cap in the polygon.  The direction is
        specified by a unit 3-vector.
        """
        return self._x

    @property
    def str(self):
        """Solid angle of this polygon (steradians).
        """
        if self._str is None:
            self._str = self.garea()
        return self._str

    def cmminf(self):
        """The index of the smallest cap in the polygon, accounting for
        negative caps and ``use_caps``.

        Returns
        -------
        :class:`int`
            Integer index of the smallest cap.
        """
        if self.ncaps == 0:
            return None
        cmmin = 2.0
        kmin = -1
        for k in range(self.ncaps):
            if is_cap_used(self.use_caps, k):
                if self.cm[k] >= 0:
                    cmk = self.cm[k]
                else:
                    cmk = 2.0 + self.cm[k]
                if cmk < cmmin:
                    cmmin = cmk
                    kmin = k
        return kmin

    def garea(self):
        """Compute the area of a polygon.

        See [1]_ for the detailed area formula, which is summarized here:

        * An empty polygon with no caps is defined to be the whole sky.
        * A polygon with one cap has area ``2*pi*self.cm``.
        * A polygon with at least one cap with an area less than :math:`2\pi`
          has an area less than :math:`2\pi`.
        * If every cap has an area greater than :math:`2\pi`, split the polygon
          into two smaller polygons and sum the two areas.

        Returns
        -------
        :class:`float`
            The area of the polygon.

        References
        ----------

        ..  [1] `Hamilton, A. J. S.; Tegmark, Max, 2004 MNRAS 349, 115
            <http://adsabs.harvard.edu/abs/2004MNRAS.349..115H>`_.
        """
        smallest_cap = self.cmminf()
        if smallest_cap is None:
            return 4.0 * np.pi
        cmmin = self.cm[smallest_cap]
        used_caps = list()
        for k in range(self.ncaps):
            if is_cap_used(self.use_caps, k):
                used_caps.append(k)
        nused_caps = len(used_caps)
        if nused_caps == 0:
            return 0.0
        if nused_caps >= 2 and cmmin > 1.0:
            npl = nused_caps + 1
        else:
            npl = nused_caps
        if npl == nused_caps:
            if nused_caps == 1:
                #
                # Short-circuit this case.
                #
                if cmmin < 0:
                    return 2.0*np.pi*(2.0+cmmin)
                else:
                    return 2.0*np.pi*cmmin
            else:
                #
                # Two or more caps, at least one has area < 2.0*pi.
                #
                return self._garea_helper()
        else:
            #
            # More than two caps, and all have area > 2.0*pi.
            # This section still needs to account for use_caps.
            #
            dpoly = self.polyn(self, smallest_cap)
            dpoly.cm[nused_caps] = cmmin / 2.0
            area1 = dpoly._garea_helper()
            dpoly.cm[nused_caps] = -1.0 * dpoly.cm[nused_caps]
            area2 = dpoly._garea_helper()
            return area1 + area2

    def _garea_helper(self):
        """Reproduces the Fortran routine garea in Mangle.

        *Placeholder for now.*

        Returns
        -------
        :class:`float`
            Area of polygon in steradians.
        """
        if self.gzeroar():
            return 0.0
        warnings.warn(("The ManglePolygon.garea() method is incomplete" +
                       "and is returning a dummy value."),
                      PydlutilsUserWarning)
        return 1.0

    def gzeroar(self):
        """If at least one cap has zero area, then the whole polygon
        has zero area.

        Returns
        -------
        :class:`bool`
            ``True`` if the area is zero.
        """
        if self.ncaps == 0:
            #
            # Catch the case of an "allsky" polygon.
            #
            return False
        return (self.cm == 0.0).any() or (self.cm <= -2.0).any()

    def copy(self):
        """Return an exact copy of the polygon.

        Returns
        -------
        :class:`ManglePolygon`
            A new polygon object.
        """
        return ManglePolygon(self)

    def add_caps(self, x, cm):
        """Create a new polygon that contains additional caps.  The caps
        are appended to the existing polygon's caps.

        Parameters
        ----------
        x : :class:`~numpy.ndarray` or :class:`~numpy.recarray`
            ``X`` values of the new caps.
        cm : :class:`~numpy.ndarray` or :class:`~numpy.recarray`
            ``CM`` values of the new caps.

        Returns
        -------
        :class:`ManglePolygon`
            A new polygon object.
        """
        ncaps = self.ncaps + cm.size
        newdata = dict()
        newdata['weight'] = self.weight
        newdata['pixel'] = self.pixel
        newdata['id'] = self.id
        newdata['use_caps'] = self.use_caps
        newdata['x'] = np.zeros((ncaps, 3), dtype=self.x.dtype)
        newdata['cm'] = np.zeros((ncaps,), dtype=self.cm.dtype)
        newdata['x'][0:self.ncaps, :] = self.x.copy()
        newdata['cm'][0:self.ncaps] = self.cm.copy()
        newdata['x'][self.ncaps:, :] = x.copy()
        newdata['cm'][self.ncaps:] = cm.copy()
        return ManglePolygon(**newdata)

    def polyn(self, other, n, complement=False):
        """Intersection of a polygon with the `n` th cap of another polygon.

        Parameters
        ----------
        other : :class:`~pydl.pydlutils.mangle.ManglePolygon`
            Polygon containing a cap to intersect the first polygon with.
        n : :class:`int`
            Index of the cap in `other`.
        complement : :class:`bool`, optional
            If ``True``, set the sign of the cm value of `other` to be
            the complement of the original value.

        Returns
        -------
        :class:`~pydl.pydlutils.mangle.ManglePolygon`
            A polygon containing the intersected caps.
        """
        x = other.x[n, :]
        sign = 1.0
        if complement:
            sign = -1.0
        cm = sign * other.cm[n]
        return self.add_caps(x, cm)


class FITS_polygon(fits.FITS_rec):
    """Handle polygons read in from a FITS file.

    This class provides aliases for the columns in a typical FITS polygons
    file.
    """
    _pkey = {'X': 'XCAPS', 'x': 'XCAPS',
             'CM': 'CMCAPS', 'cm': 'CMCAPS',
             'ID': 'IFIELD', 'id': 'IFIELD',
             'ncaps': 'NCAPS',
             'weight': 'WEIGHT',
             'pixel': 'PIXEL',
             'str': 'STR',
             'use_caps': 'USE_CAPS'}
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

    def __getitem__(self, key):
        if isinstance(key, (str,)):
            if key in self._pkey:
                return super(FITS_polygon, self).__getitem__(self._pkey[key])
        return super(FITS_polygon, self).__getitem__(key)

    def __getattr__(self, key):
        if key in self._pkey:
            return super(FITS_polygon, self).__getattribute__(self._pkey[key])
        raise AttributeError("FITS_polygon has no attribute {0}.".format(key))


def angles_to_x(points, latitude=False):
    """Convert spherical angles to unit Cartesian vectors.

    Parameters
    ----------
    points : :class:`~numpy.ndarray`
        A set of angles in the form phi, theta (in degrees).
    latitude : :class:`bool`, optional
        If ``True``, assume that the angles actually represent longitude,
        latitude, or equivalently, RA, Dec.

    Returns
    -------
    :class:`~numpy.ndarray`
        The corresponding Cartesian vectors.
    """
    npoints, ncol = points.shape
    x = np.zeros((npoints, 3), dtype=points.dtype)
    phi = np.radians(points[:, 0])
    if latitude:
        theta = np.radians(90.0 - points[:, 1])
    else:
        theta = np.radians(points[:, 1])
    st = np.sin(theta)
    x[:, 0] = np.cos(phi) * st
    x[:, 1] = np.sin(phi) * st
    x[:, 2] = np.cos(theta)
    return x


def cap_distance(x, cm, points):
    """Compute the distance from a point to a cap, and also determine
    whether the point is inside or outside the cap.

    Parameters
    ----------
    x : :class:`~numpy.ndarray` or :class:`~numpy.recarray`
        ``X`` value of the cap (3-vector).
    cm : :class:`~numpy.ndarray` or :class:`~numpy.recarray`
        ``CM`` value of the cap.
    points : :class:`~numpy.ndarray` or :class:`~numpy.recarray`
        If `points` is a 3-vector, or set of 3-vectors, then assume the point
        is a Cartesian unit vector.  If `point` is a 2-vector or set
        of 2-vectors, assume the point is RA, Dec.

    Returns
    -------
    :class:`~numpy.ndarray`
        The distance(s) to the point(s) in degrees.  If the distance is
        negative, the point is outside the cap.
    """
    npoints, ncol = points.shape
    if ncol == 2:
        xyz = angles_to_x(points, latitude=True)
    elif ncol == 3:
        xyz = points
    else:
        raise ValueError("Inappropriate shape for point!")
    dotprod = np.dot(xyz, x)
    cdist = np.degrees(np.arccos(1.0 - np.abs(cm)) - np.arccos(dotprod))
    if cm < 0:
        cdist *= -1.0
    return cdist


def circle_cap(radius, points):
    """Construct caps based on input points (directions on the unit sphere).

    Parameters
    ----------
    radius : :class:`float` or :class:`~numpy.ndarray`
        Radii of the caps in degrees.  This will become the ``CM`` value.
    points : :class:`~numpy.ndarray` or :class:`~numpy.recarray`
        If `points` is a 3-vector, or set of 3-vectors, then assume the point
        is a Cartesian unit vector.  If `point` is a 2-vector or set
        of 2-vectors, assume the point is RA, Dec.

    Returns
    -------
    :func:`tuple`
        A tuple containing ``X`` and ``CM`` values for the cap.
    """
    npoints, ncol = points.shape
    if ncol == 2:
        x = angles_to_x(points, latitude=True)
    elif ncol == 3:
        x = points.copy()
    else:
        raise ValueError("Inappropriate shape for point!")
    if isinstance(radius, float):
        radius = np.zeros((npoints,), dtype=points.dtype) + radius
    elif radius.shape == ():
        radius = np.zeros((npoints,), dtype=points.dtype) + radius
    elif radius.size == npoints:
        pass
    else:
        raise ValueError("radius does not appear to be the correct shape.")
    cm = 1.0 - np.cos(np.radians(radius))
    return (x, cm)


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


def is_in_cap(x, cm, points):
    """Are the points in a (single) given cap?

    Parameters
    ----------
    x : :class:`~numpy.ndarray` or :class:`~numpy.recarray`
        ``X`` value of the cap (3-vector).
    cm : :class:`~numpy.ndarray` or :class:`~numpy.recarray`
        ``CM`` value of the cap.
    points : :class:`~numpy.ndarray` or :class:`~numpy.recarray`
        If `points` is a 3-vector, or set of 3-vectors, then assume the point
        is a Cartesian unit vector.  If `point` is a 2-vector or set
        of 2-vectors, assume the point is RA, Dec.

    Returns
    -------
    :class:`~numpy.ndarray`
        A boolean vector giving the result for each point.
    """
    return cap_distance(x, cm, points) >= 0.0


def is_in_polygon(polygon, points, ncaps=0):
    """Are the points in a given (single) polygon?

    Parameters
    ----------
    polygon : :class:`~pydl.pydlutils.mangle.ManglePolygon`
        A polygon object.
    points : :class:`~numpy.ndarray` or :class:`~numpy.recarray`
        If `points` is a 3-vector, or set of 3-vectors, then assume the point
        is a Cartesian unit vector.  If `point` is a 2-vector or set
        of 2-vectors, assume the point is RA, Dec.
    ncaps : :class:`int`, optional
        If set, use only the first `ncaps` caps in `polygon`.

    Returns
    -------
    :class:`~numpy.ndarray`
        A boolean vector giving the result for each point.
    """
    npoints, ncol = points.shape
    p = dict()
    pmap = {'ncaps': 'NCAPS', 'use_caps': 'USE_CAPS',
            'x': 'XCAPS', 'cm': 'CMCAPS'}
    for key in pmap:
        try:
            p[key] = getattr(polygon, key)
        except AttributeError:
            p[key] = polygon[pmap[key]]
    usencaps = p['ncaps']
    if ncaps > 0:
        usencaps = min(ncaps, p['ncaps'])
    in_polygon = np.ones((npoints,), dtype=np.bool)
    for icap in range(usencaps):
        if is_cap_used(p['use_caps'], icap):
            in_polygon &= is_in_cap(p['x'][icap, :], p['cm'][icap], points)
    return in_polygon


def is_in_window(polygons, points, ncaps=0):
    """Check to see if `points` lie within a set of `polygons`.

    Parameters
    ----------
    polygons : :class:`PolygonList` or :class:`FITS_polygon`
        A set of polygons.
    points : :class:`~numpy.ndarray` or :class:`~numpy.recarray`
        If `points` is a 3-vector, or set of 3-vectors, then assume the point
        is a Cartesian unit vector.  If `point` is a 2-vector or set
        of 2-vectors, assume the point is RA, Dec.
    ncaps : :class:`int`, optional
        If set, use only the first `ncaps` caps in `polygon`.  This only
        exists to be passed to :func:`is_in_polygon`.

    Returns
    -------
    :func:`tuple`
        A tuple containing two :class:`~numpy.ndarray`.  First, a boolean
        vector giving the result for each point.  Second, an integer vector
        giving the index of the polygon that contains the point.
    """
    npoints, ncol = points.shape
    npoly = len(polygons)
    in_polygon = np.zeros((npoints, ), dtype=np.int32) - 1
    curr_polygon = 0
    while curr_polygon < npoly:
        indx_not_in = (in_polygon == -1).nonzero()[0]
        if len(indx_not_in) > 0:
            indx_in_curr_polygon = is_in_polygon(polygons[curr_polygon],
                                                 points[indx_not_in],
                                                 ncaps=ncaps)
            if indx_in_curr_polygon.any():
                in_polygon[indx_not_in[indx_in_curr_polygon]] = curr_polygon
        curr_polygon += 1
    return (in_polygon >= 0, in_polygon)


def read_fits_polygons(filename, convert=False):
    """Read a "polygon" format FITS file.

    This function returns a subclass of :class:`~astropy.io.fits.FITS_rec`
    that provides some convenient aliases for the columns of a typical
    FITS polygon file (which might be all upper-case).

    Parameters
    ----------
    filename : :class:`str`
        Name of FITS file to read.
    convert : :class:`bool`, optional
        If ``True``, convert the data to a list of :class:`ManglePolygon`
        objects.  *Caution: This could result in some data being discarded!*

    Returns
    -------
    :class:`~pydl.pydlutils.mangle.FITS_polygon` or :class:`list`
        The data contained in HDU 1 of the FITS file.
    """
    with fits.open(filename, uint=True) as hdulist:
        data = hdulist[1].data
    if convert:
        poly = PolygonList()
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
    :class:`~pydl.pydlutils.mangle.PolygonList`
        A list-like object containing
        :class:`~pydl.pydlutils.mangle.ManglePolygon` objects and
        any metadata.
    """
    with open(filename, 'rU') as ply:
        lines = ply.read().split(ply.newlines)
    try:
        npoly = int(lines[0].split()[0])
    except ValueError:
        raise PydlutilsException(("Invalid first line of {0}!  " +
                                  "Are you sure this is a Mangle " +
                                  "polygon file?").format(filename))
    p_lines = [i for i, l in enumerate(lines) if l.startswith('polygon')]
    header = lines[1:p_lines[0]]
    poly = PolygonList(header=header)
    r1 = re.compile(r'polygon\s+(\d+)\s+\(([^)]+)\):')
    mtypes = {'str': float, 'weight': float, 'pixel': int, 'caps': int}
    for p in p_lines:
        m = r1.match(lines[p])
        g = m.groups()
        pid = int(g[0])
        meta = g[1].strip().split(',')
        m1 = [m.strip().split()[1] for m in meta]
        m0 = [mtypes[m1[i]](m.strip().split()[0]) for i, m in enumerate(meta)]
        metad = dict(zip(m1, m0))
        metad['id'] = pid
        metad['x'] = list()
        metad['cm'] = list()
        for cap in lines[p+1:p+1+metad['caps']]:
            data = [float(d) for d in re.split(r'\s+', cap.strip())]
            metad['x'].append(data[0:3])
            metad['cm'].append(data[-1])
        metad['x'] = np.array(metad['x'])
        assert metad['x'].shape == (metad['caps'], 3)
        metad['cm'] = np.array(metad['cm'])
        assert metad['cm'].shape == (metad['caps'],)
        poly.append(ManglePolygon(**metad))
    return poly


def set_use_caps(polygon, index_list, add=False, tol=1.0e-10,
                 allow_doubles=False, allow_neg_doubles=False):
    """Set the bits in use_caps for a polygon.

    Parameters
    ----------
    polygon : :class:`~pydl.pydlutils.mangle.ManglePolygon`
        A polygon object.
    index_list : array-like
        A list of indices of caps to set in the polygon.  Should be no
        longer, nor contain indices greater than the number of caps
        (``polygon.ncaps``).
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
        polygon.use_caps = 0
    t2 = tol**2
    for i in index_list:
        polygon.use_caps |= (1 << index_list[i])
    if not allow_doubles:
        #
        # Check for doubles
        #
        for i in range(polygon.ncaps):
            if is_cap_used(polygon.use_caps, i):
                for j in range(i+1, polygon.ncaps):
                    if is_cap_used(polygon.use_caps, j):
                        if np.sum((polygon.x[i, :] -
                                   polygon.x[j, :])**2) < t2:
                            if ((np.absolute(polygon.cm[i] -
                                             polygon.cm[j]) < tol) or
                                ((polygon.cm[i] +
                                  polygon.cm[j]) < tol and
                                  not allow_neg_doubles)):
                                #
                                # Don't use
                                #
                                polygon.use_caps -= 1 << j
    return polygon.use_caps


def x_to_angles(points, latitude=False):
    """Convert unit Cartesian vectors to spherical angles.

    Parameters
    ----------
    points : :class:`~numpy.ndarray`
        A set of x, y, z coordinates.
    latitude : :class:`bool`, optional
        If ``True``, convert the angles to longitude,
        latitude, or equivalently, RA, Dec.

    Returns
    -------
    :class:`~numpy.ndarray`
        The corresponding spherical angles.
    """
    npoints, ncol = points.shape
    phi = np.degrees(np.arctan2(points[:, 1], points[:, 0]))
    r = (points**2).sum(1)
    theta = np.degrees(np.arccos(points[:, 2]/r))
    if latitude:
        theta = 90.0 - theta
    x = np.zeros((npoints, 2), dtype=points.dtype)
    x[:, 0] = phi
    x[:, 1] = theta
    return x


def _single_polygon(poly):
    """Given certain classes of `poly`, return a
    :class:`ManglePolygon` object.

    This function is used to allow various types of input in functions that
    take a *single* polygon as an argument.

    Parameters
    ----------
    poly : object
        An object containing data that could potentially be converted
        into a :class:`ManglePolygon`.

    Returns
    -------
    ManglePolygon
        The resulting *single* polygon object.

    Raises
    ------
    :exc:`ValueError`
        If the conversion is not possible.
    """
    if isinstance(poly, ManglePolygon):
        return poly
    if isinstance(poly, PolygonList) and len(poly) == 1:
        return poly[0]
    if isinstance(poly, FITS_polygon) and len(poly) == 1:
        return ManglePolygon(poly[0])
    if isinstance(poly, fits.fitsrec.FITS_record):
        return ManglePolygon(poly)
    raise ValueError("Can't convert input into a single polygon!")
