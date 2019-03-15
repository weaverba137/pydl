# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""This module corresponds to the sdssio directory of photoop.
"""
import os
import numpy as np


#
# Filename formats used by sdss_name and sdss_path
#
_name_formats = {
    'apObj': "{ftype}-{run:06d}-{filter}{camcol:1d}-{field:04d}.fit",
    'calibMatch': "{ftype}-{run:06d}-{camcol:1d}.fits",
    'calibPhotom': "{ftype}-{run:06d}-{camcol:1d}.fits",
    'calibPhotomGlobal': "{ftype}-{run:06d}-{camcol:1d}.fits",
    'fakeIdR': "idR-{run:06d}-{filter}{camcol:1d}-{field:04d}.fit",
    'fpAtlas': "{ftype}-{run:06d}-{camcol:1d}-{field:04d}.fit",
    'fpBIN': "{ftype}-{run:06d}-{filter}{camcol:1d}-{field:04d}.fit",
    'fpC': "{ftype}-{run:06d}-{filter}{camcol:1d}-{field:04d}.fit",
    'fpFieldStat': "{ftype}-{run:06d}-{camcol:1d}-{field:04d}.fit",
    'fpM': "{ftype}-{run:06d}-{filter}{camcol:1d}-{field:04d}.fit",
    'fpObjc': "{ftype}-{run:06d}-{camcol:1d}-{field:04d}.fit",
    'hoggObj': "{ftype}-{run:06d}-{camcol:1d}-{field:04d}.fits",
    'idFF': "{ftype}-{run:06d}-{filter}{camcol:1d}.fit",
    'idR': "{ftype}-{run:06d}-{filter}{camcol:1d}-{field:04d}.fit",
    'idRR': "{ftype}-{run:06d}-{filter}{camcol:1d}-{field:04d}.fit",
    'psBB': "{ftype}-{run:06d}-{filter}{camcol:1d}-{field:04d}.fit",
    'psFF': "{ftype}-{run:06d}-{filter}{camcol:1d}.fit",
    'psField': "{ftype}-{run:06d}-{camcol:1d}-{field:04d}.fit",
    'reObjGlobal': "{ftype}-{run:06d}-{camcol:1d}-{field:04d}.fits",
    'reObjRun': "{ftype}-{run:06d}-{camcol:1d}-{field:04d}.fits",
    'reObjTmp': "{ftype}-{run:06d}-{camcol:1d}-{field:04d}.fits",
    'tsField': "{ftype}-{run:06d}-{camcol:1d}-{rerun}-{field:04d}.fit",
    }


_path_formats = {
    'apObj': "{redux}/{rerun}/{run:d}/objcs/{camcol:1d}",
    'calibMatch': "{redux}/{rerun}/{run:d}/nfcalib",
    'calibPhotom': "{redux}/{rerun}/{run:d}/nfcalib",
    'calibPhotomGlobal': "{calib}/{rerun}/{run:d}/nfcalib",
    'fakeIdR': "{data}/{run:d}/fake_fields/{camcol:1d}",
    'fpAtlas': "{redux}/{rerun}/{run:d}/objcs/{camcol:1d}",
    'fpBIN': "{redux}/{rerun}/{run:d}/objcs/{camcol:1d}",
    'fpC': "{redux}/{rerun}/{run:d}/objcs/{camcol:1d}",
    'fpFieldStat': "{redux}/{rerun}/{run:d}/objcs/{camcol:1d}",
    'fpM': "{redux}/{rerun}/{run:d}/objcs/{camcol:1d}",
    'fpObjc': "{redux}/{rerun}/{run:d}/objcs/{camcol:1d}",
    'hoggObj': "{redux}/{rerun}/{run:d}/objcs/{camcol:1d}",
    'idFF': "{redux}/{rerun}/{run:d}/objcs/{camcol:1d}",
    'idR': "{data}/{run:d}/fields/{camcol:1d}",
    'idRR': "{data}/{run:d}/fields/{camcol:1d}",
    'psBB': "{redux}/{rerun}/{run:d}/objcs/{camcol:1d}",
    'psFF': "{redux}/{rerun}/{run:d}/objcs/{camcol:1d}",
    'psField': "{redux}/{rerun}/{run:d}/objcs/{camcol:1d}",
    'reObjGlobal': "{resolve}/{rerun}/{run:d}/resolve/{camcol:1d}",
    'reObjRun': "{redux}/{rerun}/{run:d}/resolve/{camcol:1d}",
    'reObjTmp': "{resolve}/{rerun}/{run:d}/resolve/{camcol:1d}",
    'tsField': "{redux}/{rerun}/{run:d}/calibChunks/{camcol:1d}",
    }


def filtername(f):
    """Return the name of a filter given its number.

    Parameters
    ----------
    f : :class:`int`
        The filter number.

    Returns
    -------
    :class:`str`
        The corresponding filter name.

    Examples
    --------
    >>> filtername(0)
    'u'
    """
    if isinstance(f, (str,)):
        return f
    fname = ('u', 'g', 'r', 'i', 'z')
    return fname[f]


def filternum(filt='foo'):
    """Return index number for SDSS filters either from a number or name.

    Parameters
    ----------
    filt : :class:`str`
        The filter name.

    Returns
    -------
    :class:`int`
        The corresponding filter number

    Raises
    ------
    :exc:`KeyError`
        If `filt` is not a valid filter name.

    Examples
    --------
    >>> filternum('g')
    1
    """
    if filt == 'foo':
        return list(range(5))
    else:
        filters = {'u': 0, 'g': 1, 'r': 2, 'i': 3, 'z': 4}
        return filters[filt]


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


def sdss_name(ftype, run, camcol, field, rerun='', thisfilter='r',
              no_path=False):
    """Return the name of an SDSS data file including path.

    Parameters
    ----------
    ftype : :class:`str`
        The general type of the file, for example ``'reObj'``
    run : :class:`int`
        The run number.
    camcol : :class:`int`
        The camcol number.
    field : :class:`int`
        The field number
    rerun : :class:`str`, optional
        If necessary, set the rerun name using this argument.
    thisfilter : :class:`int` or :class:`str`, optional
        If necessary, set the filter using this argument.
    no_path : :class:`bool`, optional
        Normally, sdss_name returns the full path.  If `no_path` is ``True``,
        only the basename of the file is returned.

    Returns
    -------
    :class:`str`
        The full file name, normally including the full path.

    Raises
    ------
    :exc:`KeyError`
        If the file type is unknown.
    """
    if ftype == 'reObj':
        if 'PHOTO_RESOLVE' in os.environ:
            myftype = 'reObjGlobal'
        else:
            myftype = 'reObjRun'
    else:
        myftype = ftype
    thisfilter = filtername(thisfilter)
    indict = {'ftype': myftype, 'run': run, 'camcol': camcol, 'field': field,
              'filter': thisfilter, 'rerun': rerun}
    try:
        fullname = _name_formats[myftype].format(**indict)
    except KeyError:
        raise KeyError("Unknown FTYPE = {0}".format(myftype))
    if not no_path:
        datadir = sdss_path(myftype, run, camcol, rerun)
        fullname = os.path.join(datadir, fullname)
    return fullname


def sdss_path(ftype, run, camcol=0, rerun=''):
    """Return the path name for SDSS data assuming SAS directory structure.

    Parameters
    ----------
    ftype : :class:`str`
        The general type of the file, for example ``'reObj'``
    run : :class:`int`
        The run number.
    camcol : :class:`int`, optional
        If necessary, set the camcol number using this argument.
    rerun : :class:`str`, optional
        If necessary, set the rerun name using this argument.

    Returns
    -------
    :class:`str`
        The directory in which file `ftype` lives.

    Raises
    ------
    :exc:`KeyError`
        If the file type is unknown.

    """
    indict = {
        'run': run,
        'camcol': camcol,
        'rerun': rerun,
        'calib': os.getenv('PHOTO_CALIB'),
        'data': os.getenv('PHOTO_DATA'),
        'photoobj': os.getenv('BOSS_PHOTOOBJ'),
        'redux': os.getenv('PHOTO_REDUX'),
        'resolve': os.getenv('PHOTO_RESOLVE'),
        'sky': os.getenv('PHOTO_SKY'),
        'sweep': os.getenv('PHOTO_SWEEP'),
        }
    try:
        datadir = _path_formats[ftype].format(**indict)
    except KeyError:
        raise KeyError("Unknown FTYPE = {0}".format(ftype))
    return datadir


def sdssflux2ab(flux, magnitude=False, ivar=False):
    """Convert the SDSS calibrated fluxes (magnitudes) into AB fluxes
    (magnitudes).

    Parameters
    ----------
    flux : :class:`numpy.ndarray`
        Array of calibrated fluxes or SDSS magnitudes with 5 columns,
        corresponding to the 5 filters *u*, *g*, *r*, *i*, *z*.
    magnitude : :class:`bool`, optional
        If set to ``True``, then assume `flux` are SDSS magnitudes instead of
        linear flux units.
    ivar : :class:`numpy.ndarray`, optional
        If set, the input fluxes are actually inverse variances.

    Returns
    -------
    :class:`numpy.ndarray`
        Array of fluxes or magnitudes on the AB system.

    Notes
    -----
    Uses the conversions posted by D.Hogg (sdss-calib/845)::

        u(AB,2.5m) = u(2.5m) - 0.042
        g(AB,2.5m) = g(2.5m) + 0.036
        r(AB,2.5m) = r(2.5m) + 0.015
        i(AB,2.5m) = i(2.5m) + 0.013
        z(AB,2.5m) = z(2.5m) - 0.002
    """
    #
    # Correction vector, adjust this as necessary
    #
    correction = np.array([-0.042, 0.036, 0.015, 0.013, -0.002])
    rows, cols = flux.shape
    abflux = flux.copy()
    if magnitude:
        for i in range(rows):
            abflux[i, :] += correction
    else:
        factor = 10.0**(-correction/2.5)
        if ivar:
            factor = 1.0/factor**2
        for i in range(rows):
            abflux[i, :] *= factor
    return abflux
