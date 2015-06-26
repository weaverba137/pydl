# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def fpoly(x,m):
    """Compute the first `m` simple polynomials.

    Parameters
    ----------
    x : array_like
    m : int
        The number of simple polynomials to compute.  For example, if
        ``m = 3``, x**0, x**1 and x**2 will be computed.

    Returns
    -------
    fpoly : array_like
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
