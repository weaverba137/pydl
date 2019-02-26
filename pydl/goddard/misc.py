# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""This module corresponds to the goddard/misc directory in idlutils.
"""
# Needed for Python 2 because there is a math.py module in the same directory.
from math import pi


def cirrange(ang, radians=False):
    """Convert an angle larger than 360 degrees to one less than 360 degrees.

    Parameters
    ----------
    ang : :class:`float` or array-like
        Angle to convert.  If the angle is in radians, the `radians` argument
        should be set.
    radians : class:`bool`, optional
        If ``True``, the input angle is in radians, and the output will be
        between zero and 2π.

    Returns
    -------
    :class:`float` or array-like
        Angle in the restricted range.

    Examples
    --------
    >>> from pydl.goddard.misc import cirrange
    >>> cirrange(-270.0)
    90.0
    """
    if radians:
        cnst = pi * 2.0
    else:
        cnst = 360.0
    #
    # The modulo operator automatically deals with negative values
    #
    return ang % cnst
