# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-


def spec_path(plate, path=None, topdir=None, run2d=None):
    """Return the directory containing spPlate files.

    Parameters
    ----------
    plate : :class:`int` or :class:`numpy.ndarray`
        The plate(s) to examine.
    path : :class:`str`, optional
        If set, `path` becomes the full path for every plate. In other words,
        it completely short-circuits this function.
    topdir : :class:`str`, optional
        Used to override the value of :envvar:`BOSS_SPECTRO_REDUX`.
    run2d : :class:`str`, optional
        Used to override the value of :envvar:`RUN2D`.

    Returns
    -------
    :class:`list`
        A list of directories, one for each plate.

    Raises
    ------
    KeyError
        If environment variables are not supplied.
    """
    from os import environ
    from os.path import join
    from numpy import array
    from astropy.extern.six import integer_types
    if isinstance(plate, integer_types) or plate.shape == ():
        platevec = array([plate], dtype='i4')
    else:
        platevec = plate
    if path is None:
        if topdir is None:
            topdir = environ['BOSS_SPECTRO_REDUX']
        if run2d is None:
            run2d = environ['RUN2D']
    paths = list()
    for p in platevec:
        if path is not None:
            paths.append(path)
        else:
            paths.append(join(topdir, run2d, str(p)))
    return paths
