# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def sdss_path(ftype, run, camcol=0, rerun=''):
    """Return the path name for SDSS data assuming SAS directory structure.

    Parameters
    ----------
    ftype : str
        The general type of the file, for example ``'reObj'``
    run : int
        The run number.
    camcol : int, optional
        If necessary, set the camcol number using this argument.
    rerun : str, optional
        If necessary, set the rerun name using this argument.

    Returns
    -------
    sdss_name : str
        The directory in which file `ftype` lives.

    Raises
    ------
    KeyError
        If the file type is unknown.

    """
    from os import getenv
    from . import _path_formats
    indict = {
        'run':run,
        'camcol':camcol,
        'rerun':rerun,
        'calib':getenv('PHOTO_CALIB'),
        'data':getenv('PHOTO_DATA'),
        'photoobj':getenv('BOSS_PHOTOOBJ'),
        'redux':getenv('PHOTO_REDUX'),
        'resolve':getenv('PHOTO_RESOLVE'),
        'sky':getenv('PHOTO_SKY'),
        'sweep':getenv('PHOTO_SWEEP'),
        }
    try:
        datadir = _path_formats[ftype].format(**indict)
    except KeyError:
        raise KeyError("Unknown FTYPE = {0}".format(ftype))
    return datadir
