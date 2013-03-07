# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def fchebyshev(x,m):
    """Compute a Chebyshev polynomial.

    Parameters
    ----------
    x : array_like
    m : int

    Returns
    -------
    fchebyshev : array_like
    """
    import numpy as np
    from scipy.special import chebyt
    if isinstance(x,np.ndarray):
        n = x.size
    else:
        n = 1
    if m < 1:
        raise ValueError('Order of Chebyshev polynomial must be at least 1.')
    leg = np.ones((m,n),dtype='d')
    if m >= 2:
        leg[1,:] = x
    if m >= 3:
        for k in range(2,m):
            leg[k,:] = np.polyval(chebyt(k),x)
    return leg

