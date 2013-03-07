# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def read_ds_cooling(fname,logT=None):
    """Read in Dopita & Sutherland 1993 cooling function.

    Parameters
    ----------
    fname : { 'm-00.cie', 'm-05.cie', 'm+05.cie', 'm-10.cie', 'm-15.cie', 'm-20.cie', 'm-30.cie', 'mzero.cie' }
        Name of the data file to read.
    logT : ndarray, optional
        If provided, values will be interpolated to the provided values.  If
        not provided, the values in the data files will be returned.

    Returns
    -------
    read_ds_cooling : tuple
        A tuple containing logT and loglambda, respectively.

    Raises
    ------
    URLError
        If the data file could not be retrieved.
    ValueError
        If the input file name is invalid.

    Notes
    -----
    Retrieves data from http://www.mso.anu.edu.au/~ralph/data/cool/ rather
    than storing data locally.

    Examples
    --------
    >>> pydl.pydlutils.cooling.read_ds_cooling('m-15.cie')
    """
    from urllib2 import urlopen
    from numpy import interp
    from astropy.io import ascii
    baseurl = 'http://www.mso.anu.edu.au/~ralph/data/cool/'
    if fname not in ('m-00.cie', 'm-05.cie', 'm+05.cie', 'm-10.cie', 'm-15.cie', 'm-20.cie', 'm-30.cie', 'mzero.cie'):
        raise ValueError('Invalid value for data file: {0}'.format(fname))
    with urlopen(baseurl+fname) as coolingfile:
        coolingfile = coolingfile.read()
    data = ascii.read(coolingfile.split('\n')[2:],delimiter='\t')
    if logT is None:
        return (data['log(T)'].data, data['log(lambda net)'].data)
    else:
        loglambda = interp(logT, data['log(T)'].data, data['log(lambda net)'].data)
        return (logT,loglambda)
