# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
#
__doctest_skip__ = ['read_ds_cooling']
#
def read_ds_cooling(fname,logT=None):
    """Read in Dopita & Sutherland 1993 cooling function.

    Parameters
    ----------
    fname : { 'm-00.cie', 'm-05.cie', 'm+05.cie', 'm-10.cie', 'm-15.cie', 'm-20.cie', 'm-30.cie', 'mzero.cie' }
        Name of the data file to read.
    logT : :class:`numpy.ndarray`, optional
        If provided, values will be interpolated to the provided values.  If
        not provided, the values in the data files will be returned.

    Returns
    -------
    read_ds_cooling : :func:`tuple`
        A tuple containing `logT` and `loglambda`, respectively.

    Raises
    ------
    URLError
        If the data file could not be retrieved.
    ValueError
        If the input file name is invalid.

    Notes
    -----
    Retrieves data from http://www.mso.anu.edu.au/~ralph/data/cool/ and caches
    the data locally.

    Examples
    --------
    >>> from pydl.pydlutils.cooling import read_ds_cooling
    >>> logT, loglambda = read_ds_cooling('m-15.cie')
    >>> logT[0:5]
    array([ 4.  ,  4.05,  4.1 ,  4.15,  4.2 ])
    >>> loglambda[0:5]
    array([-26.  , -24.66, -23.52, -22.62, -22.11])
    """
    from astropy.utils.data import download_file
    from numpy import interp
    from astropy.io import ascii
    baseurl = 'http://www.mso.anu.edu.au/~ralph/data/cool/'
    if fname not in ('m-00.cie', 'm-05.cie', 'm+05.cie', 'm-10.cie', 'm-15.cie', 'm-20.cie', 'm-30.cie', 'mzero.cie'):
        raise ValueError('Invalid value for data file: {0}'.format(fname))
    filename = download_file(baseurl+fname,cache=True,show_progress=False)
    with open(filename) as coolingfile:
        coolingfile = coolingfile.read()
    data = ascii.read(coolingfile.split('\n')[2:],delimiter='\t')
    if logT is None:
        return (data['log(T)'].data, data['log(lambda net)'].data)
    else:
        loglambda = interp(logT, data['log(T)'].data, data['log(lambda net)'].data)
        return (logT,loglambda)
