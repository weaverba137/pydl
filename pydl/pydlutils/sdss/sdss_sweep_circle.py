# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def sdss_sweep_circle(ra,dec,radius,stype='star',allobj=False):
    """Read the SDSS datasweep files and return objects around a location.

    Parameters
    ----------
    ra, dec : float
        The sky location to search, J2000 degrees.
    radius : float
        The radius around `ra`, `dec` to search.
    stype : str, optional
        The type of object to search, 'star', 'gal' or 'sky'.  The default is 'star'.
    allobj : bool, optional
        If set to ``True``, return all objects found, not just SURVEY_PRIMARY.

    Returns
    -------
    sdss_sweep_circle : object
        The value of the bitmask name(s).

    Notes
    -----
    Assumes that the sweep files exist in ``$PHOTO_SWEEP`` and that index files
    have been created.
    """
    import numpy as np
    from os import getenv
    from os.path import join
    from astropy.io import fits
    from ... import uniq
    from .. import PydlutilsException
    from ..spheregroup import spherematch
    from . import sweep_cache
    #
    # Check values
    #
    if stype not in ('star','gal','sky'):
        raise ValueError('Invalid type {0}!'.format(stype))
    sweepdir = getenv('PHOTO_SWEEP')
    if sweepdir is None:
        raise PydlutilsException('PHOTO_SWEEP is not set!')
    #
    # Read the index
    #
    if sweep_cache[stype] is None:
        indexfile = join(sweepdir,'datasweep-index-{0}.fits'.format(stype))
        with fits.open(indexfile) as f:
            sweep_cache[stype] = f[1].data
    index = sweep_cache[stype]
    #
    # Match
    #
    m1,m2,d12 = spherematch(np.array([ra]),np.array([dec]),index['RA'],index['DEC'],radius+0.36,maxmatch=0)
    if len(m2) == 0:
        return None
    if not allobj:
        w = index['NPRIMARY'][m2] > 0
        if w.any():
            m2 = m2[w]
        else:
            return None
    #
    # Maximum number of objects
    #
    if allobj:
        n = index['IEND'][m2] - index['ISTART'][m2] + 1
        ntot = (np.where(n>0,n,np.zeros(n.shape,dtype=n.dtype))).sum()
    else:
        ntot = index['NPRIMARY'][m2].sum()
    #
    # Find unique run + camcol
    #
    rc = index['RUN'][m2]*6 + index['CAMCOL'][m2] - 1
    isort = rc.argsort()
    iuniq = uniq(rc[isort])
    istart = 0
    return iuniq
