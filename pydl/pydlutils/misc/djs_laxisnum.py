# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def djs_laxisnum(dims,iaxis=0):
    """Returns an integer array where each element of the array is set equal to its index number in the specified axis.

    Parameters
    ----------
    dims : :class:`list`
        Dimensions of the array to return.
    iaxis : :class:`int`, optional
        Index along this dimension.

    Returns
    -------
    djs_laxisnum : :class:`numpy.ndarray`
        An array of indexes with ``dtype=int32``.

    Raises
    ------
    ValueError
        If `iaxis` is greater than or equal to the number of dimensions, or
        if number of dimensions is greater than three.

    Notes
    -----
    For two or more dimensions, there is no difference between this routine and
    :func:`~pydl.pydlutils.misc.djs_laxisgen`.

    Examples
    --------
    >>> from pydl.pydlutils.misc import djs_laxisnum
    >>> djs_laxisnum([4,4])
    array([[0, 0, 0, 0],
           [1, 1, 1, 1],
           [2, 2, 2, 2],
           [3, 3, 3, 3]], dtype=int32)
    """
    from numpy import zeros
    ndimen = len(dims)
    result = zeros(dims,dtype='i4')
    if ndimen == 1:
        pass
    elif ndimen == 2:
        if iaxis == 0:
            for k in range(dims[0]):
                result[k,:] = k
        elif iaxis == 1:
            for k in range(dims[1]):
                result[:,k] = k
        else:
            raise ValueError("Bad value for iaxis: {0:d}".format(iaxis))
    elif ndimen == 3:
        if iaxis == 0:
            for k in range(dims[0]):
                result[k,:,:] = k
        elif iaxis == 1:
            for k in range(dims[1]):
                result[:,k,:] = k
        elif iaxis == 2:
            for k in range(dims[2]):
                result[:,:,k] = k
        else:
            raise ValueError("Bad value for iaxis: {0:d}".format(iaxis))
    else:
        raise ValueError("{0:d} dimensions not supported.".format(ndimen))
    return result
