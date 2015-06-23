# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def fchebyshev(x,m):
    """Compute the first `m` Chebyshev polynomials, but modified to allow a
    split in the baseline at ``x=0``.  The intent is to allow a model fit where
    a constant term is different for positive and negative `x`.

    Parameters
    ----------
    x : array_like
    m : int
        The number of Chebyshev polynomials to compute.  For example, if
        ``m = 4``, x > 0, T_0(x), T_1(x) and T_2(x) will be computed.

    Returns
    -------
    fchebyshev : array_like
    """
    import numpy as np
    if isinstance(x,np.ndarray):
        n = x.size
    else:
        n = 1
    if m < 2:
        raise ValueError('Order of polynomial must be at least 2.')
    leg = np.ones((m,n),dtype='d')
    leg[0,:] = (x >= 0).astype(x.dtype)
    if m >= 3:
        leg[2,:] = x
        for k in range(2,m):
            leg[k,:] = leg[k-1,:]*x
            leg[k,:] = 2.0 * x * leg[k-1,:] - leg[k-2,:]
    return leg
