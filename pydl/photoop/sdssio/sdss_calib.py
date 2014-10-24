# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def sdss_calib(run,camcol,field,rerun='',**kwargs):
    """Read photometric calibration solutions from calibPhotom or calibPhotomGlobal files.

    Parameters
    ----------
    run : int
        Photo run number
    camcol : int
        Camcol number
    field : int
        Field number
    rerun : str, optional
        Photometric reduction number, as a string.

    Returns
    -------
    sdss_calib : dict
        A dictionary containing the 'NMGYPERCOUNT' keyword.

    Notes
    -----
    Currently, this is just a placeholder.
    """
    return { 'NMGYPERCOUNT':1.0 }
