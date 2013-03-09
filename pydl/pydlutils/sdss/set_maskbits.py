# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def set_maskbits(idlutils_version='v5_5_8'):
    """Populate the maskbits cache.

    Parameters
    ----------
    idlutils_version : str, optional
        Fetch the sdssMaskbits.par file corresponding to this idlutils version.

    Returns
    -------
    None

    Raises
    ------
    URLError
        If the data file could not be retrieved.
    """
    from ..yanny import yanny
    from urllib2 import urlopen
    from astropy.utils.data import download_file
    if idlutils_version == 'trunk' or idlutils_version.startswith('branches/'):
        iversion = idlutils_version
    else:
        iversion = 'tags/'+idlutils_version
    baseurl = 'http://www.sdss3.org/svn/repo/idlutils/{0}/data/sdss/sdssMaskbits.par'.format(iversion)
    filename = download_file(baseurl,cache=True)
    #par = urlopen(baseurl)
    maskfile = yanny(filename)
    #par.close()
    #
    # Parse the file & cache the results in maskbits
    #
    maskbits = dict()
    for k in range(maskfile.size('MASKBITS')):
        if maskfile['MASKBITS']['flag'][k] in maskbits:
            maskbits[maskfile['MASKBITS']['flag'][k]][maskfile['MASKBITS']['label'][k]] = maskfile['MASKBITS']['bit'][k]
        else:
            maskbits[maskfile['MASKBITS']['flag'][k]] = {maskfile['MASKBITS']['label'][k]:maskfile['MASKBITS']['bit'][k]}
    return maskbits


