# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def fchebyshev_split(x,m):
    """Compute the first `m` Chebyshev polynomials, but modified to allow a
    split in the baseline at :math:`x=0`.  The intent is to allow a model fit where
    a constant term is different for positive and negative `x`.

    Parameters
    ----------
    x : array-like
        Compute the Chebyshev polynomials at these abscissa values.
    m : :class:`int`
        The number of Chebyshev polynomials to compute.  For example, if
        :math:`m = 3`, :math:`T_0 (x)`, :math:`T_1 (x)` and :math:`T_2 (x)` will be computed.

    Returns
    -------
    fchebyshev_split : :class:`numpy.ndarray`
    """
    import numpy as np
    if isinstance(x,np.ndarray):
        n = x.size
    else:
        n = 1
    if m < 2:
        raise ValueError('Order of polynomial must be at least 2.')
    try:
        dt = x.dtype
    except AttributeError:
        dt = np.float64
    leg = np.ones((m,n),dtype=dt)
    try:
        leg[0,:] = (x >= 0).astype(x.dtype)
    except AttributeError:
        leg[0,:] = np.double(x >= 0)
    if m > 2:
        leg[2,:] = x
    if m > 3:
        for k in range(3,m):
            leg[k,:] = 2.0 * x * leg[k-1,:] - leg[k-2,:]
    return leg
