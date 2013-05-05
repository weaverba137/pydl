# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def sdss_calibv():
    """Return calibration for velocities from pix/frame to deg/day.

    Parameters
    ----------
    None

    Returns
    -------
    sdss_calibv : float
        The conversion from pixels per frame to degrees per day

    Notes
    -----
    Assumes frame time difference of 71.72 seconds and pixel scale of 0.396 arcsec, both fixed.
    Also note that observations of the same part of sky from adjacent bands are separated by *two* frame numbers,
    so we multiply by a factor two.
    """
    pixscale = 0.396  # arcsec
    ftime = 71.72    # seconds
    pixframe2degday= 2.0*pixscale/(3600.0) * (3600.0)*24.0/ftime
    return pixframe2degday
