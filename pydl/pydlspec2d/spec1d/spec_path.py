# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from __future__ import print_function
from astropy.extern import six


def spec_path(plate,path=None,topdir=None,run2d=None):
    """Return the path to spPlate files.

    Parameters
    ----------
    plate : int or ndarray
        The plate(s) to examine.
    path : str, optional
        If set, `path` becomes the full path for every plate. In other words,
        it completely short-circuits this function.
    topdir : str, optional
        Used to override the value of BOSS_SPECTRO_REDUX.
    run2d : str, optional
        Used to override the value of RUN2D.

    Returns
    -------
    spec_path : list
        A list of paths, one for each plate.
    """
    from os import getenv
    from os.path import join
    from numpy import array
    if isinstance(plate, six.integer_types) or plate.shape == ():
        platevec = array([plate],dtype='i4')
    else:
        platevec = plate
    if path is None:
        if topdir is None:
            topdir = getenv('BOSS_SPECTRO_REDUX')
        if run2d is None:
            run2d = getenv('RUN2D')
    paths = list()
    for p in platevec:
        if path is not None:
            paths.append(path)
        else:
            paths.append(join(topdir,run2d,str(p)))
    return paths

