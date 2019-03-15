# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""This module corresponds to the window directory in photoop.
"""
import os
import time
from warnings import warn
import numpy as np
from astropy import log
from astropy.io import fits
from . import PhotoopException
from .sdssio import sdss_name, sdss_calib
from ..pydlutils.mangle import set_use_caps
from ..pydlutils.sdss import sdss_flagval


def sdss_score(flist, silent=True):
    """Score a list of imaging fields from zero to one.

    Parameters
    ----------
    flist : :class:`~astropy.io.fits.HDUList`
        Opened FITS file, typically ``window_flist.fits``.
    silent : :class:`bool`, optional
        If ``False``, print extra information.

    Returns
    -------
    :class:`numpy.ndarray`
        A vector of scores, one for each row of the FITS file.
    """
    lat = 32.780361
    lng = 360.0 - 105.820417
    tzone = 7
    scores = 1
    #
    # Read the PHOTO status bits
    #
    if not silent:
        log.info('Setting PHOTO status bits')
    t1 = time.time()
    nlist = flist[1].header.get('NAXIS2')
    fdata = flist[1].data
    for k in range(nlist):
        if not silent and ((k % 1000) == 0):
            log.info("Setting PHOTO status {0:d} {1:d}".format(k, nlist))
        thisfile = sdss_name('fpFieldStat', fdata.field('RUN')[k],
                             fdata.field('CAMCOL')[k], fdata.field('FIELD')[k],
                             fdata.field('RERUN')[k])
        try:
            fpfield = fits.open(thisfile)
        except IOError:
            warn("Retrying fpFieldStat file for RUN = {0:d} CAMCOL = {1:d} FIELD = {2:d} RERUN = {3}".format(fdata.field('RUN')[k],
                 fdata.field('CAMCOL')[k], fdata.field('FIELD')[k],
                 fdata.field('RERUN')[k]))
            try:
                fpfield = fits.open(thisfile)
            except IOError:
                warn("Bad fpFieldStat file for RUN = {0:d} CAMCOL = {1:d} FIELD = {2:d} RERUN = {3}".format(fdata.field('RUN')[k],
                     fdata.field('CAMCOL')[k], fdata.field('FIELD')[k],
                     fdata.field('RERUN')[k]))
                fdata.field('PHOTO_STATUS')[k] = -1
                if not silent:
                    log.info('Trying tsField instead.')
                thisfile = sdss_name('tsField', fdata.field('RUN')[k],
                                     fdata.field('CAMCOL')[k],
                                     fdata.field('FIELD')[k],
                                     fdata.field('RERUN')[k])
                try:
                    tsfield = fits.open(thisfile)
                except IOError:
                    warn('Bad tsField file.')
                else:
                    if not silent:
                        log.info('tsField found, using frames_status.')
                    fdata.field('PHOTO_STATUS')[k] = tsfield[1].data.field('frames_status')[0]
            else:
                fdata.field('PHOTO_STATUS')[k] = fpfield[1].data.field('status')[0]
        else:
            fdata.field('PHOTO_STATUS')[k] = fpfield[1].data.field('status')[0]
    if not silent:
        log.info("Time to set PHOTO status = {0:f} sec".format(time.time()-t1))
    #
    # Read in the PSP status
    #
    if not silent:
        log.info('Setting PSP status bits')
    t2 = time.time()
    for k in range(nlist):
        if not silent and ((k % 1000) == 0):
            log.info("Setting PSP status {0:d} {1:d}".format(k, nlist))
        thisfile = sdss_name('psField', fdata.field('RUN')[k],
                             fdata.field('CAMCOL')[k], fdata.field('FIELD')[k],
                             fdata.field('RERUN')[k])
        try:
            psfield = fits.open(thisfile)
        except IOError:
            warn("Retrying psField file for RUN = {0:d} CAMCOL = {1:d} FIELD = {2:d} RERUN = {3}".format(fdata.field('RUN')[k],
                 fdata.field('CAMCOL')[k], fdata.field('FIELD')[k],
                 fdata.field('RERUN')[k]))
            try:
                psfield = fits.open(thisfile)
            except IOError:
                warn("Bad psField file for RUN = {0:d} CAMCOL = {1:d} FIELD = {2:d} RERUN = {3}".format(fdata.field('RUN')[k],
                     fdata.field('CAMCOL')[k], fdata.field('FIELD')[k],
                     fdata.field('RERUN')[k]))
                fdata.field('PSP_STATUS')[k] = -1
                fdata.field('PSF_FWHM')[k] = -1
                fdata.field('SKYFLUX')[k] = -1
        pixscale = 0.396 * np.sqrt(fdata.field('XBIN')[k]**2 +
                                   fdata.field('YBIN')[k]**2)/np.sqrt(2.0)
        calibinfo = sdss_calib(fdata.field('RUN')[k],
                               fdata.field('CAMCOL')[k],
                               fdata.field('FIELD')[k],
                               fdata.field('RERUN')[k], **kwargs)
        fdata.field('PSP_STATUS')[k] = psfield[6].data.field('status')[0]
        fdata.field('PSF_FWHM')[k] = psfield[6].data.field('psf_width')[0]
        fdata.field('SKYFLUX')[k] = (psfield[6].data.field('sky')[0] *
                                     calibinfo['NMGYPERCOUNT']/pixscale**2)
    if not silent:
        log.info("Time to set PSP status = {0:f} sec".format(time.time()-t2))
    #
    # Decide if each field exists in all 5 bands.
    #
    bad_bits = sdss_flagval('image_status', ['bad_rotator', 'bad_astrom',
                            'bad_focus', 'shutters'])
    if 'ignoreframesstatus' in kwargs:
        ignoreframesstatus = np.zeros(fdata.field('PHOTO_STATUS').shape) == 0
    else:
        ignoreframesstatus = np.zeros(fdata.field('PHOTO_STATUS').shape) == 1
    qexist = (fdata.field('PHOTO_STATUS') == 0) | ignoreframesstatus
    for k in range(5):
        qexist &= (fdata.field('IMAGE_STATUS')[:, k] & bad_bits) == 0
    #
    # Decide if each field is phtometric in all 5 bands.
    #
    unphot_bits = sdss_flagval('image_status', ['cloudy', 'unknown',
                               'ff_petals', 'dead_ccd', 'noisy_ccd'])
    qphot = fdata.field('SUN_ANGLE') < -12
    for k in range(5):
        qphot &= (fdata.field('IMAGE_STATUS')[:, k] & unphot_bits) == 0
    for k in range(5):
        qphot &= (((fdata.field('PSP_STATUS')[:, k] & 31) <= 2) |
                  (fdata.field('XBIN') > 1) | ignoreframesstatus)
    #
    # Now set the score for each field
    #
    sensitivity = (0.7 / (fdata.field('PSF_FWHM')[:, 2] *
                   np.sqrt(fdata.field('SKYFLUX')[:, 2]))) < 0.4
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
        warn("Changing NaN scores for {0:d} fields to zero.".format(len(ibad)))
        fdata.field('SCORE')[ibad] = 0
    return fdata.field('SCORE')


def window_read(flist=False, rescore=False, blist=False, bcaps=False,
                balkans=False, findx=False, bindx=False):
    """Read window files in :envvar:`PHOTO_RESOLVE`.

    Parameters
    ----------
    flist : :class:`bool`, optional
        If ``True``, read the ``window_flist.fits`` file.
    rescore : :class:`bool`, optional
        If `flist` is ``True``, look for ``window_flist_rescore.fits``, and run
        :func:`~pydl.photoop.window.window_score` if it is not found.
    blist : :class:`bool`, optional
        If ``True``, read the ``window_blist.fits`` file.
    bcaps : :class:`bool`, optional
        If ``True``, read the ``window_bcaps.fits`` file.
    balkans : :class:`bool`, optional
        If ``True``, construct the balkans from the ``window_blist.fits`` and
        ``window_bcaps.fits`` files.
    findx : :class:`bool`, optional
        If ``True``, read the ``window_findx.fits`` file.
    bindx : :class:`bool`, optional
        If ``True``, read the ``window_bindx.fits`` file.

    Returns
    -------
    :class:`dict`
        A dictionary containing the requested window data.
    """
    try:
        resolve_dir = os.environ['PHOTO_RESOLVE']
    except KeyError:
        raise PhotoopException(('You have not set the environment variable ' +
                                'PHOTO_RESOLVE!'))
    r = dict()
    if flist:
        flist_file = os.path.join(resolve_dir, 'window_flist.fits')
        if rescore:
            flist_file = os.path.join(resolve_dir, 'window_flist_rescore.fits')
            if not os.path.exists(flist_file):
                window_score(rescore=rescore)
        with fits.open(flist_file) as fit:
            r['flist'] = fit[1].data
    if blist or balkans:
        blist_file = os.path.join(resolve_dir, 'window_blist.fits')
        with fits.open(blist_file) as fit:
            r['blist'] = fit[1].data
    if bcaps or balkans:
        bcaps_file = os.path.join(resolve_dir, 'window_bcaps.fits')
        with fits.open(bcaps_file) as fit:
            r['bcaps'] = fit[1].data
    if findx:
        findx_file = os.path.join(resolve_dir, 'window_findx.fits')
        with fits.open(findx_file) as fit:
            r['findx'] = fit[1].data
    if bindx:
        bindx_file = os.path.join(resolve_dir, 'window_bindx.fits')
        with fits.open(bindx_file) as fit:
            r['bindx'] = fit[1].data
    if balkans:
        #
        # Copy blist data to balkans
        #
        r['balkans'] = r['blist'].copy()
        r['balkans']['caps'] = {'X': list(), 'CM': list()}
        r['balkans']['use_caps'] = np.zeros(r['balkans']['ICAP'].shape,
                                            dtype=np.uint64)
        if not blist:
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
        if not bcaps:
            del r['bcaps']
    return r


def window_score(rescore=False):
    """For uber-resolve, score all the fields from zero to one.

    Parameters
    ----------
    rescore : :class:`bool`, optional
        If `rescore` is ``True``, then write a new file 'window_flist_rescore.fits'
        rather than over-writing the file 'window_flist.fits'
    """
    #
    # Be certain not to use global calibrations
    #
    try:
        calib_dir_save = os.environ['PHOTO_CALIB']
    except KeyError:
        raise PhotoopException(('You have not set the environment variable ' +
                                'PHOTO_CALIB!'))
    del os.environ['PHOTO_CALIB']
    #
    # Read the file
    #
    try:
        resolve_dir = os.environ['PHOTO_RESOLVE']
    except KeyError:
        raise PhotoopException(('You have not set the environment variable ' +
                                'PHOTO_RESOLVE!'))
    filename = os.path.join(resolve_dir, 'window_flist.fits')
    if rescore:
        fitsmode = 'readonly'
    else:
        fitsmode = 'update'
    try:
        flist = fits.open(filename, mode=fitsmode)
    except OSError:
        raise PhotoopException('Unable to read FLIST file.')
    #
    # Construct the scores filling in the values to FLIST.SCORE
    #
    flist[1].data['SCORE'][:] = sdss_score(flist)
    if rescore:
        flist.writeto(os.path.join(resolve_dir, 'window_flist_rescore.fits'))
    flist.close()
    #
    # Restore the PHOTO_CALIB variable
    #
    os.environ['PHOTO_CALIB'] = calib_dir_save
    return
