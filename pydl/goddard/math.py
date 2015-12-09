# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""This module corresponds to the goddard/math directory in idlutils.
"""
from __future__ import absolute_import


def flegendre(x, m):
    """Compute the first `m` Legendre polynomials.

    Parameters
    ----------
    x : array-like
        Compute the Legendre polynomials at these abscissa values.
    m : :class:`int`
        The number of Legendre polynomials to compute.  For example, if
        :math:`m = 3`, :math:`P_0 (x)`, :math:`P_1 (x)` and :math:`P_2 (x)`
        will be computed.

    Returns
    -------
    :class:`numpy.ndarray`
    """
    import numpy as np
    from scipy.special import legendre
    if isinstance(x, np.ndarray):
        n = x.size
    else:
        n = 1
    if m < 1:
        raise ValueError('Number of Legendre polynomials must be at least 1.')
    try:
        dt = x.dtype
    except AttributeError:
        dt = np.float64
    leg = np.ones((m, n), dtype=dt)
    if m >= 2:
        leg[1, :] = x
    if m >= 3:
        for k in range(2, m):
            leg[k, :] = np.polyval(legendre(k), x)
    return leg
