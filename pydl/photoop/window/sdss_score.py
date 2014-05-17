# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def sdss_score(flist,**kwargs):
    """Score a list of imaging fields from zero to one.
    """
    import time
    import numpy as np
    import astropy.io.fits as pyfits
    from ..sdssio import sdss_name, sdss_calib
    from ...pydlutils.sdss import sdss_flagval
    if 'silent' not in kwargs:
        kwargs['silent'] = True
    lat = 32.780361
    lng = 360.0 - 105.820417
    tzone = 7
    scores = 1
    #
    # Read the PHOTO status bits
    #
    if not kwargs['silent']:
        print('Setting PHOTO status bits')
    t1 = time.time()
    nlist = flist[1].header.get('NAXIS2')
    fdata = flist[1].data
    for k in range(nlist):
        if not kwargs['silent'] and ((k % 1000) == 0):
            print("Setting PHOTO status {0:d} {1:d}".format(k, nlist))
        thisfile = sdss_name('fpFieldStat',fdata.field('RUN')[k],
            fdata.field('CAMCOL')[k],fdata.field('FIELD')[k],
            fdata.field('RERUN')[k])
        try:
            fpfield = pyfits.open(thisfile)
        except IOError:
            print("Warning: Retrying fpFieldStat file for RUN = {0:d} CAMCOL = {1:d} FIELD = {2:d} RERUN = {3}".format(fdata.field('RUN')[k],
                fdata.field('CAMCOL')[k],fdata.field('FIELD')[k],
                fdata.field('RERUN')[k]))
            try:
                fpfield = pyfits.open(thisfile)
            except IOError:
                print("Warning: Bad fpFieldStat file for RUN = {0:d} CAMCOL = {1:d} FIELD = {2:d} RERUN = {3}".format(fdata.field('RUN')[k],
                    fdata.field('CAMCOL')[k],fdata.field('FIELD')[k],
                    fdata.field('RERUN')[k]))
                fdata.field('PHOTO_STATUS')[k] = -1
                if not kwargs['silent']:
                    print('Trying tsField instead.')
                thisfile = sdss_name('tsField', fdata.field('RUN')[k],
                    fdata.field('CAMCOL')[k],fdata.field('FIELD')[k],
                    fdata.field('RERUN')[k])
                try:
                    tsfield = pyfits.open(thisfile)
                except IOError:
                    print('Warning: Bad tsField file.')
                else:
                    if not kwargs['silent']:
                        print('tsField found, using frames_status.')
                    fdata.field('PHOTO_STATUS')[k] = tsfield[1].data.field('frames_status')[0]
            else:
                fdata.field('PHOTO_STATUS')[k] = fpfield[1].data.field('status')[0]
        else:
            fdata.field('PHOTO_STATUS')[k] = fpfield[1].data.field('status')[0]
    if not kwargs['silent']:
        print("Time to set PHOTO status = {0:f} sec".format(time.time()-t1))
    #
    # Read in the PSP status
    #
    if not kwargs['silent']:
        print('Setting PSP status bits')
    t2 = time.time()
    for k in range(nlist):
        if not kwargs['silent'] and ((k % 1000) == 0):
            print("Setting PSP status {0:d} {1:d}".format(k, nlist))
        thisfile = sdss_name('psField',fdata.field('RUN')[k],
            fdata.field('CAMCOL')[k],fdata.field('FIELD')[k],
            fdata.field('RERUN')[k])
        try:
            psfield = pyfits.open(thisfile)
        except IOError:
            print("Warning: Retrying psField file for RUN = {0:d} CAMCOL = {1:d} FIELD = {2:d} RERUN = {3}".format(fdata.field('RUN')[k],
                fdata.field('CAMCOL')[k],fdata.field('FIELD')[k],
                fdata.field('RERUN')[k]))
            try:
                psfield = pyfits.open(thisfile)
            except IOError:
                print("Warning: Bad psField file for RUN = {0:d} CAMCOL = {1:d} FIELD = {2:d} RERUN = {3}".format(fdata.field('RUN')[k],
                    fdata.field('CAMCOL')[k],fdata.field('FIELD')[k],
                    fdata.field('RERUN')[k]))
                fdata.field('PSP_STATUS')[k] = -1
                fdata.field('PSF_FWHM')[k] = -1
                fdata.field('SKYFLUX')[k] = -1
        pixscale = 0.396 * np.sqrt(fdata.field('XBIN')[k]**2 + fdata.field('YBIN')[k]**2)/np.sqrt(2.0)
        calibinfo = sdss_calib(fdata.field('RUN')[k],
            fdata.field('CAMCOL')[k],fdata.field('FIELD')[k],
            fdata.field('RERUN')[k],**kwargs)
        fdata.field('PSP_STATUS')[k] = psfield[6].data.field('status')[0]
        fdata.field('PSF_FWHM')[k] = psfield[6].data.field('psf_width')[0]
        fdata.field('SKYFLUX')[k] = psfield[6].data.field('sky')[0] * calibinfo['NMGYPERCOUNT']/ pixscale**2
    if not kwargs['silent']:
        print("Time to set PSP status = {0:f} sec".format(time.time()-t2))
    #
    # Decide if each field exists in all 5 bands.
    #
    bad_bits = sdss_flagval('image_status',
        ['bad_rotator','bad_astrom','bad_focus','shutters'])
    if 'ignoreframesstatus' in kwargs:
        ignoreframesstatus = np.zeros(fdata.field('PHOTO_STATUS').shape) == 0
    else:
        ignoreframesstatus = np.zeros(fdata.field('PHOTO_STATUS').shape) == 1
    qexist = (fdata.field('PHOTO_STATUS') == 0) | ignoreframesstatus
    for k in range(5):
        qexist &= (fdata.field('IMAGE_STATUS')[:,k] & bad_bits) == 0
    #
    # Decide if each field is phtometric in all 5 bands.
    #
    unphot_bits = sdss_flagval('image_status',
        ['cloudy','unknown','ff_petals','dead_ccd','noisy_ccd'])
    qphot = fdata.field('SUN_ANGLE') < -12
    for k in range(5):
        qphot &= (fdata.field('IMAGE_STATUS')[:,k] & unphot_bits) == 0
    for k in range(5):
        qphot &= (((fdata.field('PSP_STATUS')[:,k] & 31) <= 2) |
            (fdata.field('XBIN') > 1) | ignoreframesstatus)
    #
    # Now set the score for each field
    #
    sensitivity = (0.7 / (fdata.field('PSF_FWHM')[:,2] * np.sqrt(fdata.field('SKYFLUX')[:,2]))) < 0.4
    fdata.field('SCORE')[:] = qexist * (0.1 + 0.5*qphot + sensitivity)
    ibinned = np.find(fdata.field('XBIN') > 1)
    if len(ibinned) > 0:
        fdata.field('SCORE')[ibinned] *= 0.1
    #
    # Look for any NaN values, which could happen for example if there
    # is a corrupted psField file and PSF_FWHM or SKYFLUX was negative.
    #
    ibad = np.find(~np.isfinite(fdata.field('SCORE')))
    if len(ibad) > 0:
        print("Warning: Changing NaN scores for {0:d} fields to zero.".format(len(ibad)))
        fdata.field('SCORE')[ibad] = 0
    return fdata.field('SCORE')
