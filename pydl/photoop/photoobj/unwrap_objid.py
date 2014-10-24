# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def unwrap_objid(objid):
    """Unwrap CAS-style objID into run, camcol, field, id, rerun.

    See pydl.pydlutils.sdss.sdss_objid() for details on how the bits
    within an objID are assigned.

    Parameters
    ----------
    objid : numpy.ndarray
        An array containing 64-bit integers or strings.  If strings are passed,
        they will be converted to integers internally.

    Returns
    -------
    unwrap_objid : numpy.recarray
        A record array with the same length as objid, with the columns
        'run', 'camcol', 'frame', 'id', 'rerun', 'skyversion'.

    Notes
    -----
    For historical reasons, the inverse of this function,
    pydl.pydlutils.sdss.sdss_objid() is not in the same namespace as this
    function.

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
    from numpy import bitwise_and, int64, recarray, string_, unicode_
    if objid.dtype.type is string_ or objid.dtype.type is unicode_:
        tempobjid = objid.astype(int64)
    elif objid.dtype.type is int64:
        tempobjid = objid.copy()
    else:
        raise ValueError('Unrecognized type for objid!')
    unwrap = recarray(objid.shape,dtype=[('skyversion','i4'),('rerun','i4'),('run','i4'),
        ('camcol','i4'),('frame','i4'),('id','i4')])
    unwrap.skyversion = bitwise_and(tempobjid >> 59, 2**4 - 1)
    unwrap.rerun = bitwise_and(tempobjid >> 48, 2**11 - 1)
    unwrap.run = bitwise_and(tempobjid >> 32, 2**16 - 1)
    unwrap.camcol = bitwise_and(tempobjid >> 29, 2**3 - 1)
    # unwrap.firstfield = bitwise_and(tempobjid >> 28, 2**1 - 1)
    unwrap.frame = bitwise_and(tempobjid >> 16, 2**12 - 1)
    unwrap.id = bitwise_and(tempobjid, 2**16 - 1)
    return unwrap
