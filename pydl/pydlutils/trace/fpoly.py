# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def fpoly(x,m):
    """Compute the first `m` simple polynomials.

    Parameters
    ----------
    x : array-like
        Compute the simple polynomials at these abscissa values.
    m : :class:`int`
        The number of simple polynomials to compute.  For example, if
        :math:`m = 3`, :math:`x^0`, :math:`x^1` and :math:`x^2` will be computed.

    Returns
    -------
    fpoly : :class:`numpy.ndarray`
    """
    import numpy as np
    if isinstance(x,np.ndarray):
        n = x.size
    else:
        n = 1
    if m < 1:
        raise ValueError('Order of polynomial must be at least 1.')
    try:
        dt = x.dtype
    except AttributeError:
        dt = np.float64
    leg = np.ones((m,n),dtype=dt)
    if m >= 2:
        leg[1,:] = x
    if m >= 3:
        for k in range(2,m):
            leg[k,:] = leg[k-1,:]*x
    return leg
