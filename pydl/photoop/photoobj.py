# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""This module corresponds to the photoobj directory in photoop.
"""


def sdss_calibv():
    """Return calibration for velocities from pix/frame to deg/day.

    Returns
    -------
    :class:`float`
        The conversion from pixels per frame to degrees per day

    Notes
    -----
    Assumes frame time difference of 71.72 seconds and pixel scale of
    0.396 arcsec, both fixed. Also note that observations of the same part of
    sky from adjacent bands are separated by *two* frame numbers,
    so we multiply by a factor two.
    """
    pixscale = 0.396   # arcsec
    ftime = 71.72    # seconds
    pixframe2degday = 2.0*pixscale/(3600.0) * (3600.0)*24.0/ftime
    return pixframe2degday


def unwrap_objid(objid):
    """Unwrap CAS-style objID into run, camcol, field, id, rerun.

    See :func:`~pydl.pydlutils.sdss.sdss_objid` for details on how the bits
    within an objID are assigned.

    Parameters
    ----------
    objid : :class:`numpy.ndarray`
        An array containing 64-bit integers or strings.  If strings are passed,
        they will be converted to integers internally.

    Returns
    -------
    :class:`numpy.recarray`
        A record array with the same length as objid, with the columns
        'run', 'camcol', 'frame', 'id', 'rerun', 'skyversion'.

    Notes
    -----
    For historical reasons, the inverse of this function,
    :func:`~pydl.pydlutils.sdss.sdss_objid` is not in the same namespace as
    this function.

    'frame' is used instead of 'field' because record arrays have a method
    of the same name.

    Examples
    --------
    >>> from numpy import array
    >>> from pydl.photoop.photoobj import unwrap_objid
    >>> unwrap_objid(array([1237661382772195474]))
    rec.array([(2, 301, 3704, 3, 91, 146)],
          dtype=[('skyversion', '<i4'), ('rerun', '<i4'), ('run', '<i4'), ('camcol', '<i4'), ('frame', '<i4'), ('id', '<i4')])
    """
    import numpy as np
    if objid.dtype.type is np.string_ or objid.dtype.type is np.unicode_:
        tempobjid = objid.astype(np.int64)
    elif objid.dtype.type is np.int64:
        tempobjid = objid.copy()
    else:
        raise ValueError('Unrecognized type for objid!')
    unwrap = np.recarray(objid.shape,
                         dtype=[('skyversion', 'i4'), ('rerun', 'i4'),
                                ('run', 'i4'), ('camcol', 'i4'),
                                ('frame', 'i4'), ('id', 'i4')])
    unwrap.skyversion = np.bitwise_and(tempobjid >> 59, 2**4 - 1)
    unwrap.rerun = np.bitwise_and(tempobjid >> 48, 2**11 - 1)
    unwrap.run = np.bitwise_and(tempobjid >> 32, 2**16 - 1)
    unwrap.camcol = np.bitwise_and(tempobjid >> 29, 2**3 - 1)
    # unwrap.firstfield = np.bitwise_and(tempobjid >> 28, 2**1 - 1)
    unwrap.frame = np.bitwise_and(tempobjid >> 16, 2**12 - 1)
    unwrap.id = np.bitwise_and(tempobjid, 2**16 - 1)
    return unwrap
