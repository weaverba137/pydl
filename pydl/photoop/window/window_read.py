# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def window_read(**kwargs):
    """Read window files in $PHOTO_RESOLVE.
    """
    from os import getenv
    from os.path import join
    from . import window_score
    from .. import PhotoopException
    from ...goddard.fits import mrdfits
    from ...pydlutils.mangle import set_use_caps
    import numpy as np
    resolve_dir = getenv('PHOTO_RESOLVE')
    if resolve_dir is None:
        raise PhotoopException('You have not set the environment variable PHOTO_RESOLVE!')
    if 'silent' not in kwargs:
        kwargs['silent'] = True
    r = dict()
    if 'flist' in kwargs:
        if 'rescore' in kwargs:
            flist,h = mrdfits(join(resolve_dir,'window_flist_rescore.fits'),1,**kwargs)
            if flist is None:
                #
                # This will be called if window_flist_rescore.fits doesn't exist.
                #
                window_score()
                flist,h = mrdfits(join(resolve_dir,'window_flist_rescore.fits'),1,**kwargs)
        else:
            flist,h = mrdfits(join(resolve_dir,'window_flist.fits'),1,**kwargs)
        r['flist'] = flist
    if 'blist' in kwargs or 'balkans' in kwargs:
        r['blist'],h = mrdfits(join(resolve_dir,'window_blist.fits'),1,**kwargs)
    if 'bcaps' in kwargs or 'balkans' in kwargs:
        r['bcaps'],h = mrdfits(join(resolve_dir,'window_bcaps.fits'),1,**kwargs)
    if 'findx' in kwargs:
        r['findx'],h = mrdfits(join(resolve_dir,'window_findx.fits'),1,**kwargs)
    if 'bindx' in kwargs:
        r['bindx'],h = mrdfits(join(resolve_dir,'window_bindx.fits'),1,**kwargs)
    if 'balkans' in kwargs:
        #
        # Copy blist data to balkans
        #
        r['balkans'] = r['blist'].copy()
        r['balkans']['caps'] = { 'X':list(), 'CM':list() }
        r['balkans']['use_caps'] = pylab.zeros(r['balkans']['ICAP'].shape,dtype=pylab.uint64)
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
