# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from __future__ import print_function
#
#
#
def get_juldate(seconds=None):
    """Returns the current Julian date.

    Uses the MJD trick & adds the offset to get JD.

    Parameters
    ----------
    seconds : int or float, optional
        Time in seconds since the UNIX epoch.  This should only be used
        for testing.

    Returns
    -------
    get_juldate : float
        The Julian Day number as a floating point number.

    Notes
    -----
    Do not use this function if high precision is required, or if you are
    concerned about the distinction between UTC & TAI.
    """
    from time import time
    if seconds is None:
        t = time()
    else:
        t = seconds
    mjd = t/86400.0 + 40587.0
    return mjd + 2400000.5
#
#
#
def get_juldate_main(): # pragma: no cover
    """Allow this module to be run in scripts.
    """
    jd = get_juldate()
    print(jd)
    return 0
