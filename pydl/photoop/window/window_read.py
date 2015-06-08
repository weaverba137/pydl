# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def window_read(**kwargs):
    """Read window files in $PHOTO_RESOLVE.
    """
    from os import getenv
    from os.path import exists, join
    from . import window_score
    from .. import PhotoopException
    from ...pydlutils.mangle import set_use_caps
    from astropy.io import fits
    import numpy as np
    resolve_dir = getenv('PHOTO_RESOLVE')
    if resolve_dir is None:
        raise PhotoopException('You have not set the environment variable PHOTO_RESOLVE!')
    if 'silent' not in kwargs:
        kwargs['silent'] = True
    r = dict()
    if 'flist' in kwargs:
        if 'rescore' in kwargs:
            flist_file = join(resolve_dir,'window_flist_rescore.fits')
            if not exists(rescore_file):
                #
                # This will be called if window_flist_rescore.fits doesn't exist.
                #
                window_score()
        else:
            flist_file = join(resolve_dir,'window_flist.fits')
        with fits.open(rescore_file) as fit:
            r['flist'] = fit[1].data
    if 'blist' in kwargs or 'balkans' in kwargs:
        blist_file = join(resolve_dir,'window_blist.fits')
        with fits.open(balkan_file) as fit:
            r['blist'] = fit[1].data
    if 'bcaps' in kwargs or 'balkans' in kwargs:
        bcaps_file = join(resolve_dir,'window_bcaps.fits')
        with fits.open(bcaps_file) as fit:
            r['bcaps'] = fit[1].data
    if 'findx' in kwargs:
        findx_file = join(resolve_dir,'window_findx.fits')
        with fits.open(findx_file) as fit:
            r['findx'] = fit[1].data
    if 'bindx' in kwargs:
        bindx_file = join(resolve_dir,'window_bindx.fits')
        with fits.open(bindx_file) as fit:
            r['bindx'] = fit[1].data
    if 'balkans' in kwargs:
        #
        # Copy blist data to balkans
        #
        r['balkans'] = r['blist'].copy()
        r['balkans']['caps'] = { 'X':list(), 'CM':list() }
        r['balkans']['use_caps'] = np.zeros(r['balkans']['ICAP'].shape,dtype=np.uint64)
        if 'blist' not in kwargs:
            del r['blist']
        #
        # Copy bcaps data into balkans
        #
        for k in range(len(r['balkans']['ICAP'])):
            r['balkans']['caps']['X'].append(r['bcaps']['X'][balkans['ICAP'][k]:balkans['ICAP'][k]+balkans['NCAPS'][k]])
            r['balkans']['caps']['CM'].append(r['bcaps']['CM'][balkans['ICAP'][k]:balkans['ICAP'][k]+balkans['NCAPS'][k]])
            r['balkans']['use_caps'][k] = set_use_caps(
                r['balkans']['caps']['X'][k],
                r['balkans']['caps']['CM'][k],
                r['balkans']['use_caps'][k],
                allow_doubles=True)
        if 'bcaps' not in kwargs:
            del r['bcaps']
    return r
