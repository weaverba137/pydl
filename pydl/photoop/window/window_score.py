# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def window_score(**kwargs):
    """
    For uber-resolve, score all the fields from zero to one.  If
    'rescore' is set, then write a new file 'window_flist_rescore.fits'
    rather than over-writing the file 'window_flist.fits'
    """
    # import time
    from os import environ
    from astropy.io import fits as pyfits
    from . import sdss_score
    from .. import PhotoopException
    # t0 = time.time()
    #
    # Be certain not to use global calibrations
    #
    try:
        calib_dir_save = environ['PHOTO_CALIB']
    except KeyError:
        raise PhotoopException('You have not set the environment variable PHOTO_CALIB!')
    del environ['PHOTO_CALIB']
    #
    # Read the file
    #
    try:
        resolve_dir = environ['PHOTO_RESOLVE']
    except KeyError:
        raise PhotoopException('You have not set the environment variable PHOTO_RESOLVE!')
    filename = join(resolve_dir,'window_flist.fits')
    if 'rescore' in kwargs:
        fitsmode = 'readonly'
    else:
        fitsmode = 'update'
    try:
        flist = pyfits.open(filename,mode=fitsmode)
    except IOError:
        raise PhotoopException('Unable to read FLIST file.')
    #
    # Construct the scores filling in the values to FLIST.SCORE
    #
    flist.field('SCORE')[:] = sdss_score(flist)
    if 'rescore' in kwargs:
        flist.writeto(join(resolve_dir,'window_flist_rescore.fits'))
    flist.close()
    #
    # Restore the PHOTO_CALIB variable
    #
    environ['PHOTO_CALIB'] = calib_dir_save
    # print "Elapsed time = %f sec" % (time.time()-t0,)
    return
