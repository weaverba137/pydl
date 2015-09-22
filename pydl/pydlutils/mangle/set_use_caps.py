# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from astropy.extern import six

if six.PY3:
    # This works, but wouldn't `int` intead of `long` be OK to use here?
    # http://python3porting.com/differences.html#long
    long = int


def set_use_caps(x,cm,polygon_use_caps,add=False,tol=1.0e-10,allow_doubles=False,allow_neg_doubles=False):
    """Set the bits in use_caps for a polygon.

    Parameters
    ----------
    x : array-like
        Polygon `x` value.
    cm : array-like
        Polygon `cm` value.
    polygon_use_caps : :class:`int`
        Input value of use_caps.
    add : :class:`bool`, optional
        If ``True``, don't initialize the use_caps value to zero, use the value of
        `polygon_use_caps` instead.
    tol : :class:`float`, optional
        Tolerance used to determine whether two caps are identical.
    allow_doubles : :class:`bool`, optional
        Normally, this routine automatically sets use_caps such that no
        two caps with use_caps set are identical.
    allow_neg_doubles : :class:`bool`, optional
        Normally, two caps that are identical except for the sign of `cm`
        would be set unused.  This inhibits that behaviour.

    Returns
    -------
    set_use_caps : :class:`long`
        Value of use_caps.
    """
    import numpy as np
    from . import is_cap_used
    if add:
        use_caps = long(polygon_use_caps)
    else:
        use_caps = long(0)
    t2 = tol**2
    use_caps |= long(2)**len(cm) - long(1)
    if not allow_doubles:
        #
        # Check for doubles
        #
        for i in range(len(cm)):
            if is_cap_used(use_caps,i):
                for j in range(i+1,len(cm)):
                    if is_cap_used(use_caps,j):
                        if np.sum((x[i]-x[j])**2) < t2:
                            if ((np.absolute(cm[i]-cm[j]) < tol) or
                                ((cm[i] + cm[j]) < tol and not allow_neg_doubles)):
                                #
                                # Don't use
                                #
                                use_caps -= long(2)**j
    return use_caps
