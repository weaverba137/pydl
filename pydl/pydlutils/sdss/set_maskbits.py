# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def set_maskbits(idlutils_version='v5_5_8',maskbits_file=None):
    """Populate the maskbits cache.

    Parameters
    ----------
    idlutils_version : :class:`str`, optional
        Fetch the sdssMaskbits.par file corresponding to this idlutils version.
    maskbits_file : :class:`str`, optional
        Use an explicit file instead of downloading the official version.
        This should only be used for tests.

    Returns
    -------
    set_maskbits : :class:`dict`
        A dictionary of bitmasks suitable for caching.

    Raises
    ------
    URLError
        If the data file could not be retrieved.
    """
    from ..yanny import yanny
    from astropy.utils.data import download_file
    if maskbits_file is None: # pragma: no cover
        if idlutils_version == 'trunk' or idlutils_version.startswith('branches/'):
            iversion = idlutils_version
        else:
            iversion = 'tags/'+idlutils_version
        baseurl = 'http://www.sdss3.org/svn/repo/idlutils/{0}/data/sdss/sdssMaskbits.par'.format(iversion)
        filename = download_file(baseurl,cache=True,show_progress=False)
    else:
        filename = maskbits_file
    maskfile = yanny(filename)
    #
    # Parse the file & cache the results in maskbits
    #
    maskbits = dict()
    for k in range(maskfile.size('MASKBITS')):
        if maskfile['MASKBITS']['flag'][k] in maskbits:
            maskbits[maskfile['MASKBITS']['flag'][k]][maskfile['MASKBITS']['label'][k]] = maskfile['MASKBITS']['bit'][k]
        else:
            maskbits[maskfile['MASKBITS']['flag'][k]] = {maskfile['MASKBITS']['label'][k]:maskfile['MASKBITS']['bit'][k]}
    if 'MASKALIAS' in maskfile:
        for k in range(maskfile.size('MASKALIAS')):
            maskbits[maskfile['MASKALIAS']['alias'][k]] = maskbits[maskfile['MASKALIAS']['flag'][k]].copy()
    return maskbits
