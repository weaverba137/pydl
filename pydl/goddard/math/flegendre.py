# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def flegendre(x,m):
    """Compute the first `m` Legendre polynomials.

    Parameters
    ----------
    x : array_like
    m : int
        The number of Legendre polynomials to compute.  For example, if
        ``m = 3``, P_0(x), P_1(x) and P_2(x) will be computed.

    Returns
    -------
    flegendre : array_like
    """
    import numpy as np
    from scipy.special import legendre
    if isinstance(x,np.ndarray):
        n = x.size
    else:
        n = 1
    if m < 1:
        raise ValueError('Number of Legendre polynomials must be at least 1.')
    try:
        dt = x.dtype
    except AttributeError:
        dt = np.float64
    leg = np.ones((m,n),dtype=dt)
    if m >= 2:
        leg[1,:] = x
    if m >= 3:
        for k in range(2,m):
            leg[k,:] = np.polyval(legendre(k),x)
    return leg
