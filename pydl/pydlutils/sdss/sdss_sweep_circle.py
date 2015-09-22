# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def sdss_sweep_circle(ra,dec,radius,stype='star',allobj=False):
    """Read the SDSS datasweep files and return objects around a location.

    Parameters
    ----------
    ra, dec : :class:`float`
        The sky location to search, J2000 degrees.
    radius : :class:`float`
        The radius around `ra`, `dec` to search.
    stype : :class:`str`, optional
        The type of object to search, 'star', 'gal' or 'sky'.  The default is 'star'.
    allobj : :class:`bool`, optional
        If set to ``True``, return all objects found, not just SURVEY_PRIMARY.

    Returns
    -------
    sdss_sweep_circle : `numpy.ndarray`
        The data extracted from the sweep files.

    Notes
    -----
    Assumes that the sweep files exist in :envvar:`PHOTO_SWEEP` and that index files
    have been created.
    """
    import numpy as np
    from os import getenv
    from os.path import join
    from astropy.io import fits
    from ... import uniq
    from .. import PydlutilsException
    from ..spheregroup import spherematch
    from . import sweep_cache, sdss_flagval
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
    ira = np.array([ra])
    idec = np.array([dec])
    m1,m2,d12 = spherematch(index['RA'],index['DEC'],ira,idec,radius+0.36,maxmatch=0)
    if len(m1) == 0:
        return None
    if not allobj:
        w = index['NPRIMARY'][m1] > 0
        if w.any():
            m1 = m1[w]
        else:
            return None
    #
    # Maximum number of objects
    #
    if allobj:
        n = index['IEND'][m1] - index['ISTART'][m1] + 1
        ntot = (np.where(n>0,n,np.zeros(n.shape,dtype=n.dtype))).sum()
    else:
        ntot = index['NPRIMARY'][m1].sum()
    #
    # Find unique run + camcol
    #
    rc = index['RUN'][m1]*6 + index['CAMCOL'][m1] - 1
    isort = rc.argsort()
    iuniq = uniq(rc[isort])
    istart = 0
    objs = None
    nobjs = 0
    for i in range(len(iuniq)):
        iend=iuniq[i]
        icurr=isort[istart:iend]
        #
        # Determine which file and range of rows
        #
        run = index['RUN'][m1[icurr[0]]]
        camcol = index['CAMCOL'][m1[icurr[0]]]
        rerun = index['RERUN'][m1[icurr[0]]]
        fields = index[m1[icurr]]
        ist = fields['ISTART'].min()
        ind = fields['IEND'].max()
        if ind >= ist:
            #
            # Read in the rows of that file
            #
            swfile = join(getenv('PHOTO_SWEEP'),rerun,
                'calibObj-{0:06d}-{1:1d}-{2}.fits.gz'.format(int(run),int(camcol),stype))
            with fits.open(swfile) as f:
                tmp_objs = f[1].data[ist:ind]
            if tmp_objs.size > 0:
                #
                # Keep only objects within the desired radius
                #
                tm1,tm2,d12 = spherematch(tmp_objs['RA'],tmp_objs['DEC'],ira,idec,radius,maxmatch=0)
                if len(tm1) > 0:
                    tmp_objs=tmp_objs[tm1]
                    #
                    # Keep only SURVEY_PRIMARY objects by default
                    #
                    if not allobj:
                        w = (tmp_objs['RESOLVE_STATUS'] & sdss_flagval('RESOLVE_STATUS','SURVEY_PRIMARY')) > 0
                        if w.any():
                            tmp_objs = tmp_objs[w]
                        else:
                            tmp_objs = None
                    if tmp_objs is not None:
                        if objs is None:
                            objs = np.zeros(ntot,dtype=tmp_objs.dtype)
                        objs[nobjs:nobjs+tmp_objs.size] = tmp_objs
                        nobjs += tmp_objs.size
        istart=iend+1
    if nobjs > 0:
        return objs[0:nobjs]
    else:
        return None
