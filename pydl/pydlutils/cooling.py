# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""This module corresponds to the cooling directory in idlutils.
"""


def read_ds_cooling(fname, logT=None):
    """Read in the `Sutherland & Dopita (1993)
    <http://adsabs.harvard.edu/abs/1993ApJS...88..253S>`_ cooling function.

    Parameters
    ----------
    fname : { 'm-00.cie', 'm-05.cie', 'm+05.cie', 'm-10.cie', 'm-15.cie', 'm-20.cie', 'm-30.cie', 'mzero.cie' }
        Name of the data file to read.
    logT : :class:`numpy.ndarray`, optional
        If provided, values will be interpolated to the provided values.  If
        not provided, the values in the data files will be returned.

    Returns
    -------
    :func:`tuple`
        A tuple containing `logT` and `loglambda`, respectively.

    Raises
    ------
    ValueError
        If the input file name is invalid.

    Notes
    -----
    The data have been retrieved from
    http://www.mso.anu.edu.au/~ralph/data/cool/ and stored in the package.

    Examples
    --------
    >>> from pydl.pydlutils.cooling import read_ds_cooling
    >>> logT, loglambda = read_ds_cooling('m-15.cie')
    >>> logT[0:5]
    array([ 4.  ,  4.05,  4.1 ,  4.15,  4.2 ])
    >>> loglambda[0:5]
    array([-26.  , -24.66, -23.52, -22.62, -22.11])
    """
    from astropy.utils.data import get_pkg_data_contents
    from numpy import interp
    from astropy.io import ascii
    if fname not in ('m-00.cie', 'm-05.cie', 'm+05.cie', 'm-10.cie',
                     'm-15.cie', 'm-20.cie', 'm-30.cie', 'mzero.cie'):
        raise ValueError('Invalid value for data file: {0}'.format(fname))
    coolingfile = get_pkg_data_contents('data/cooling/' + fname)
    data = ascii.read(coolingfile.split('\n')[2:], delimiter='\t')
    if logT is None:
        return (data['log(T)'].data, data['log(lambda net)'].data)
    else:
        loglambda = interp(logT, data['log(T)'].data,
                           data['log(lambda net)'].data)
        return (logT, loglambda)
