# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def cirrange(ang,radians=False):
    """Convert an angle larger than 360 degrees to one less than 360 degrees.

    Parameters
    ----------
    ang : float or array_like
        Angle to convert.  If the angle is in radians, the `radians` argument should be set.
    radians : bool
        If true, the input angle is in radians, and the output will be between
        zero and 2*pi.

    Returns
    -------
    cirrange : float or array_like
        Angle in the restricted range.

    Examples
    --------
    >>> from pydl.goddard.misc import cirrange
    >>> cirrange(-270.0)
    90.0
    """
    from math import pi
    if radians:
        cnst = pi * 2.0
    else:
        cnst = 360.0
    #
    # The modulo operator automatically deals with negative values
    #
    return ang % cnst
