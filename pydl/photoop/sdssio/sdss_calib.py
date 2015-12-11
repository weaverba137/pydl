# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-


def sdss_calib(run, camcol, field, rerun='', **kwargs):
    """Read photometric calibration solutions from calibPhotom or
    calibPhotomGlobal files.

    Parameters
    ----------
    run : :class:`int`
        Photo run number
    camcol : :class:`int`
        Camcol number
    field : :class:`int`
        Field number
    rerun : :class:`str`, optional
        Photometric reduction number, as a string.

    Returns
    -------
    :class:`dict`
        A dictionary containing the 'NMGYPERCOUNT' keyword.

    Notes
    -----
    Currently, this is just a placeholder.
    """
    return {'NMGYPERCOUNT': 1.0}
