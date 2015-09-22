# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def djs_laxisgen(dims,iaxis=0):
    """Returns an integer array where each element of the array is set equal to its index number along the specified axis.

    Parameters
    ----------
    dims : :class:`list`
        Dimensions of the array to return.
    iaxis : :class:`int`, optional
        Index along this dimension.

    Returns
    -------
    djs_laxisgen : :class:`numpy.ndarray`
        An array of indexes with ``dtype=int32``.

    Raises
    ------
    ValueError
        If `iaxis` is greater than or equal to the number of dimensions.

    Notes
    -----
    For two or more dimensions, there is no difference between this routine and
    :func:`~pydl.pydlutils.misc.djs_laxisnum`.

    Examples
    --------
    >>> from pydl.pydlutils.misc import djs_laxisgen
    >>> djs_laxisgen([4,4])
    array([[0, 0, 0, 0],
           [1, 1, 1, 1],
           [2, 2, 2, 2],
           [3, 3, 3, 3]], dtype=int32)
    """
    from . import djs_laxisnum
    from numpy import arange
    ndimen = len(dims)
    if ndimen == 1:
        return arange(dims[0],dtype='i4')
    return djs_laxisnum(dims,iaxis)
