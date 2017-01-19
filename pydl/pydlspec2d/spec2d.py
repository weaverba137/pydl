# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""This module corresponds to the spec2d directory in idlspec2d.
"""


def aesthetics(flux, invvar, method='traditional'):
    """Add nice values to a spectrum where it is masked.

    Parameters
    ----------
    flux : :class:`numpy.ndarray`
        The spectrum to clean up.
    invvar : :class:`numpy.ndarray`
        Inverse variance of the spectrum.
    method : { 'traditional', 'noconst', 'mean', 'damp', 'nothing' }, optional
        Apply this method to clean up the spectrum.  Default is 'traditional'.

    Returns
    -------
    :class:`numpy.ndarray`
        A cleaned-up spectrum.
    """
    import numpy as np
    from scipy.special import erf
    from ..pydlutils.image import djs_maskinterp
    from . import Pydlspec2dException
    badpts = invvar == 0
    if badpts.any():
        if method == 'traditional':
            newflux = djs_maskinterp(flux, invvar == 0, const=True)
        elif method == 'noconst':
            newflux = djs_maskinterp(flux, invvar == 0)
        elif method == 'mean':
            newflux = flux.copy()
            goodpts = invvar > 0
            newflux[~goodpts] = newflux[goodpts].mean()
        elif method == 'damp':
            l = 250  # damping length in pixels
            goodpts = invvar.nonzero()[0]
            nflux = flux.size
            mingood = goodpts.min()
            maxgood = goodpts.max()
            newflux = djs_maskinterp(flux, invvar == 0, const=True)
            pixels = np.arange(nflux, dtype='f')
            if mingood > 0:
                damp1 = float(min(mingood, l))
                newflux *= 0.5*(1.0+erf((pixels-mingood)/damp1))
            if maxgood < (nflux - 1):
                damp2 = float(min(maxgood, l))
                newflux *= 0.5*(1.0+erf((maxgood-pixels)/damp2))
        elif method == 'nothing':
            newflux = flux.copy()
        else:
            raise Pydlspec2dException("Unknown method: {0}".format(method))
        return newflux
    else:
        return flux


def combine1fiber(inloglam, objflux, newloglam, objivar=None, verbose=False,
                  **kwargs):
    """Combine several spectra of the same object, or resample a single spectrum.

    Parameters
    ----------
    inloglam : :class:`numpy.ndarray`
        Vector of log wavelength.
    objflux : :class:`numpy.ndarray`
        Input flux.
    newloglam : :class:`numpy.ndarray`
        Output wavelength pixels, vector of log wavelength.
    objivar : :class:`numpy.ndarray`, optional
        Inverse variance of the flux.
    verbose : :class:`bool`, optional
        If ``True``, set log level to DEBUG.

    Returns
    -------
    :func:`tuple` of :class:`numpy.ndarray`
        The resulting flux and inverse variance.

    Raises
    ------
    ValueError
        If input dimensions don't match.
    """
    import numpy as np
    from astropy import log
    from warnings import warn
    from . import Pydlspec2dUserWarning
    from .. import smooth
    from ..pydlutils.bspline import iterfit
    from ..pydlutils.math import djs_median
    from ..pydlutils.sdss import sdss_flagval
    #
    # Log
    #
    # log.enable_warnings_logging()
    if verbose:
        log.setLevel('DEBUG')
    #
    # Check that dimensions of inputs are valid.
    #
    npix = inloglam.size
    nfinalpix = len(newloglam)
    if objflux.shape != inloglam.shape:
        raise ValueError('Dimensions of inloglam and objflux do not agree.')
    if objivar is not None:
        if objivar.shape != inloglam.shape:
            raise ValueError('Dimensions of inloglam and objivar do not agree.')
    if 'finalmask' in kwargs:
        if kwargs['finalmask'].shape != inloglam.shape:
            raise ValueError('Dimensions of inloglam and finalmask do not agree.')
    if 'indisp' in kwargs:
        if kwargs['indisp'].shape != inloglam.shape:
            raise ValueError('Dimensions of inloglam and indisp do not agree.')
    #
    # Set defaults
    #
    EPS = np.finfo(np.float32).eps
    if 'binsz' in kwargs:
        binsz = kwargs['binsz']
    else:
        if inloglam.ndim == 2:
            binsz = inloglam[0, 1] - inloglam[0, 0]
        else:
            binsz = inloglam[1] - inloglam[0]
    if 'nord' in kwargs:
        nord = kwargs['nord']
    else:
        nord = 3
    if 'bkptbin' in kwargs:
        bkptbin = kwargs['bkptbin']
    else:
        bkptbin = 1.2 * binsz
    if 'maxsep' in kwargs:
        maxsep = kwargs['maxsep']
    else:
        maxsep = 2.0 * binsz
    if inloglam.ndim == 1:
        #
        # Set specnum = 0 for all elements
        #
        nspec = 1
        specnum = np.zeros(inloglam.shape, dtype=inloglam.dtype)
    else:
        nspec, ncol = inloglam.shape
        specnum = np.tile(np.arange(nspec), ncol).reshape(ncol, nspec).transpose()
    #
    # Use fullcombmask for modifying the pixel masks in the original input files.
    #
    fullcombmask = np.zeros(npix)
    newflux = np.zeros(nfinalpix, dtype=inloglam.dtype)
    newmask = np.zeros(nfinalpix, dtype='i4')
    newivar = np.zeros(nfinalpix, dtype=inloglam.dtype)
    newdisp = np.zeros(nfinalpix, dtype=inloglam.dtype)
    newsky = np.zeros(nfinalpix, dtype=inloglam.dtype)
    newdispweight = np.zeros(nfinalpix, dtype=inloglam.dtype)
    if objivar is None:
        nonzero = np.arange(npix, dtype='i4')
        ngood = npix
    else:
        nonzero = (objivar.ravel() > 0).nonzero()[0]
        ngood = nonzero.size
    #
    # ormask is needed to create andmask
    #
    andmask = np.zeros(nfinalpix, dtype='i4')
    ormask = np.zeros(nfinalpix, dtype='i4')
    if ngood == 0:
        #
        # In this case of no good points, set the nodata bit everywhere.
        # Also if noplug is set in the first input bit-mask, assume it
        # should be set everywhere in the output bit masks.  No other bits
        # are set.
        #
        warn('No good points!', Pydlspec2dUserWarning)
        bitval = sdss_flagval('SPPIXMASK', 'NODATA')
        if 'finalmask' in kwargs:
            bitval |= (sdss_flagval('SPPIXMASK', 'NOPLUG') *
                       (finalmask[0] & sdss_flagval('SPPIXMASK', 'NODATA')))
        andmask = andmask | bitval
        ormask = ormask | bitval
        return (newflux, newivar)
    else:
        #
        # Now let's break sorted wavelengths into groups where pixel
        # separations are larger than maxsep.
        #
        inloglam_r = inloglam.ravel()
        isort = nonzero[inloglam_r[nonzero].argsort()]
        wavesort = inloglam_r[isort]
        padwave = np.insert(wavesort, 0, wavesort.min() - 2.0*maxsep)
        padwave = np.append(padwave, wavesort.max() + 2.0*maxsep)
        ig1 = ((padwave[1:ngood+1]-padwave[0:ngood]) > maxsep).nonzero()[0]
        ig2 = ((padwave[2:ngood+2]-padwave[1:ngood+1]) > maxsep).nonzero()[0]
        if ig1.size != ig2.size:
            raise ValueError('Grouping tricks did not work!')
        #
        # Avoid flux-dependent bias when combining multiple spectra.
        # This call to djs_median contains a width that is both floating-point
        # and even, which is very strange.
        #
        if objivar is not None and objivar.ndim > 1:
            saved_objivar = objivar
            for spec in range(nspec):
                igood = (objivar[spec, :] > 0).nonzero()[0]
                if igood.size > 0:
                    # objivar[spec, igood] = djs_median(saved_objivar[spec, igood], width=100.)
                    objivar[spec, igood] = djs_median(saved_objivar[spec, igood], width=101)
        else:
            saved_objivar = None
        for igrp in range(ig1.size):
            ss = isort[ig1[igrp]:ig2[igrp]+1]
            if ss.size > 2:
                if objivar is None:
                    #
                    # Fit without variance
                    #
                    sset, bmask = iterfit(inloglam_r[ss],
                                          objflux.ravel()[ss],
                                          nord=nord, groupbadpix=True,
                                          requiren=1, bkspace=bkptbin,
                                          silent=True)
                else:
                    #
                    # Fit with variance
                    #
                    sset, bmask = iterfit(inloglam_r[ss],
                                          objflux.ravel()[ss],
                                          invvar=objivar.ravel()[ss],
                                          nord=nord, groupbadpix=True,
                                          requiren=1, bkspace=bkptbin,
                                          silent=True)
                if np.sum(np.absolute(sset.coeff)) == 0:
                    sset = None
                    bmask = np.zeros(len(ss))
                    warn('All B-spline coefficients have been set to zero!',
                         Pydlspec2dUserWarning)
            else:
                bmask = np.zeros(len(ss))
                sset = None
                warn('Not enough data for B-spline fit!', Pydlspec2dUserWarning)
            inside = ((newloglam >= (inloglam_r[ss]).min()-EPS) &
                      (newloglam <= (inloglam_r[ss]).max()+EPS)).nonzero()[0]
            #
            # It is possible for numinside to be zero, if the input data points
            # span an extremely small wavelength range, within which there are
            # no output wavelengths.
            #
            if sset is not None and len(inside) > 0:
                newflux[inside], bvalumask = sset.value(newloglam[inside])
                if bvalumask.any():
                    newmask[inside[bvalumask]] = 1
                log.debug('Masked {0:d} of {1:d} pixels.'.format((1-bmask).sum(), bmask.size))
                #
                # Determine which pixels should be masked based upon the spline
                # fit. Set the combinerej bit.
                #
                ireplace = ~bmask
                if ireplace.any():
                    #
                    # The following would replace the original flux values of
                    # masked pixels with b-spline evaluations.
                    #
                    # objflux[ss[ireplace]] = sset.value(inloglam[ss[ireplace]])
                    #
                    # Set the inverse variance of these pixels to zero.
                    #
                    if objivar is not None:
                        objivar.ravel()[ss[ireplace]] = 0.0
                        log.debug('Replaced {0:d} pixels in objivar.'.format(len(ss[ireplace])))
                    if 'finalmask' in kwargs:
                        finalmask[ss[ireplace]] = (finalmask[ss[ireplace]] |
                                                   sdss_flagval('SPPIXMASK',
                                                   'COMBINEREJ'))
            fullcombmask[ss] = bmask
        #
        # Restore objivar
        #
        if saved_objivar is not None:
            objivar = saved_objivar * (objivar > 0)
        #
        # Combine inverse variance and pixel masks.
        #
        # Start with all bits set in andmask
        #
        andmask[:] = -1
        for j in range(int(specnum.max())+1):
            these = (specnum.ravel() == j).nonzero()[0]
            if these.any():
                inbetween = ((newloglam >= inloglam_r[these].min()) &
                             (newloglam <= inloglam_r[these].max()))
                if inbetween.any():
                    jnbetween = inbetween.nonzero()[0]
                    #
                    # Conserve inverse variance by doing a linear interpolation
                    # on that quantity.
                    #
                    result = np.interp(newloglam[jnbetween], inloglam_r[these],
                                       (objivar.ravel()[these] *
                                        fullcombmask[these]))
                    #
                    # Grow the fullcombmask below to reject any new sampling
                    # containing even a partial masked pixel.
                    #
                    smask = np.interp(newloglam[jnbetween], inloglam_r[these],
                                      fullcombmask[these].astype(inloglam.dtype))
                    result *= smask >= (1.0 - EPS)
                    newivar[jnbetween] += result*newmask[jnbetween]
                lowside = np.floor((inloglam_r[these]-newloglam[0])/binsz).astype('i4')
                highside = lowside + 1
                if 'finalmask' in kwargs:
                    andmask[lowside] &= finalmask[these]
                    andmask[highside] &= finalmask[these]
                    ormask[lowside] |= finalmask[these]
                    ormask[highside] |= finalmask[these]
                #
                # Combine the dispersions + skies in the dumbest way possible
                # [sic].
                #
                if 'indisp' in kwargs:
                    newdispweight[jnbetween] += result
                    newdisp[jnbetween] += (result *
                                           np.interp(newloglam[jnbetween],
                                                     inloglam_r[these],
                                                     indisp.ravel()[these]))
                    newsky[jnbetween] += (result *
                                          np.interp(newloglam[jnbetween],
                                                    inloglam_r[these],
                                                    skyflux.ravel()[these]))
        if 'indisp' in kwargs:
            newdisp /= newdispweight + (newdispweight == 0)
            newsky /= newdispweight + (newdispweight == 0)
    #
    # Grow regions where 3 or more pixels are rejected together ???
    #
    foo = smooth(newivar, 3)
    badregion = np.absolute(foo) < EPS
    # badregion = foo == 0.0
    if badregion.any():
        warn('Growing bad pixel region, {0:d} pixels found.'.format(badregion.sum()),
             Pydlspec2dUserWarning)
        ibad = badregion.nonzero()[0]
        lowerregion = np.where(ibad-2 < 0, 0, ibad-2)
        upperregion = np.where(ibad+2 > nfinalpix-1, nfinalpix-1, ibad+2)
        newivar[lowerregion] = 0.0
        newivar[upperregion] = 0.0
    #
    # Replace NaNs in combined spectra; this should really never happen.
    #
    inff = ((~np.isfinite(newflux)) | (~np.isfinite(newivar)))
    if inff.any():
        warn('{0:d} NaNs in combined spectra.'.format(inff.sum()),
             Pydlspec2dUserWarning)
        newflux[inff] = 0.0
        newivar[inff] = 0.0
    #
    # Interpolate over masked pixels, just for aesthetic purposoes.
    #
    goodpts = newivar > 0
    if 'aesthetics' in kwargs:
        amethod = kwargs['aesthetics']
    else:
        amethod = 'traditional'
    newflux = aesthetics(newflux, newivar, method=amethod)
    # if 'interpolate' in kwargs:
    #     newflux = pydlutils.image.djs_maskinterp(newflux,~goodpts,const=True)
    # else:
    #     newflux[~goodpts] = newflux[goodpts].mean()
    if goodpts.any():
        minglam = newloglam[goodpts].min()
        maxglam = newloglam[goodpts].max()
        ibad = ((newloglam < minglam) | (newloglam > maxglam))
        if ibad.any():
            ormask[ibad] |= sdss_flagval('SPPIXMASK', 'NODATA')
            andmask[ibad] |= sdss_flagval('SPPIXMASK', 'NODATA')
    #
    # Replace values of -1 in the andmask with 0.
    #
    andmask *= (andmask != -1)
    return (newflux, newivar)


def filter_thru(flux, waveimg=None, wset=None, mask=None,
                filter_prefix='sdss_jun2001', toair=False):
    """Compute throughput in SDSS filters.

    Parameters
    ----------
    flux : array-like
        Spectral flux.
    waveimg : array-like, optional
        Full wavelength solution, with the same shape as `flux`.
    wset : :class:`~pydl.pydlutils.trace.TraceSet`, optional
        A trace set containing the wavelength solution.  Must be specified
        if `waveimg` is not specified.
    mask : array-like, optional
        Interpolate over pixels where `mask` is non-zero.
    filter_prefix : :class:`str`, optional
        Specifies a set of filter curves.
    toair : :class:`bool`, optional
        If ``True``, convert the wavelengths to air from vacuum before computing.

    Returns
    -------
    array-like
        Integrated flux in the filter bands.

    Raises
    ------
    ValueError
        If neither `waveimg` nor `wset` are set.
    """
    import numpy as np
    from astropy.io import ascii
    from ..goddard.astro import vactoair
    from ..pydlutils.image import djs_maskinterp
    from ..pydlutils.trace import traceset2xy, xy2traceset
    nTrace, nx = flux.shape
    if filter_prefix != 'sdss_jun2001':
        raise ValueError("Filters other than {0} are not available!".format('sdss_jun2001'))
    ffiles = [_get_pkg_filename_compat('data/filters/{0}_{1}_atm.dat'.format(filter_prefix, f),
                                    package='pydl.pydlutils') for f in 'ugriz']
    if waveimg is None and wset is None:
        raise ValueError("Either waveimg or wset must be specified!")
    if waveimg is None:
        pixnorm, logwave = traceset2xy(wset)
        waveimg = 10**logwave
    if toair:
        newwaveimg = vactoair(waveimg)
    else:
        newwaveimg = waveimg
    logwave = np.log10(newwaveimg)
    diffx = np.outer(np.ones((nTrace,), dtype=flux.dtype), np.arange(nx-1, dtype=flux.dtype))
    diffy = logwave[:, 1:] - logwave[:, 0:nx-1]
    diffset = xy2traceset(diffx, diffy, ncoeff=4, xmin=0, xmax=nx-1)
    pixnorm, logdiff = traceset2xy(diffset)
    logdiff = np.absolute(logdiff)
    if mask is not None:
        flux_interp = djs_maskinterp(flux, mask, iaxis=0)
    res = np.zeros((nTrace, len(ffiles)), dtype=flux.dtype)
    for i, f in enumerate(ffiles):
        filter_data = ascii.read(f, comment='#.*', names=('lam', 'respt',
                                 'resbig', 'resnoa', 'xatm'))
        filtimg = logdiff * np.interp(newwaveimg.flatten(),
                                      filter_data['lam'].data,
                                      filter_data['respt'].data).reshape(logdiff.shape)
        if mask is not None:
            res[:, i] = (flux_interp * filtimg).sum(1)
        else:
            res[:, i] = (flux * filtimg).sum(1)
        sumfilt = filtimg.sum(1)
        res[:, i] = res[:, i] / (sumfilt + (sumfilt <= 0).astype(sumfilt.dtype))
    return res


def _get_pkg_filename_compat(filename, package):
    """Astropy 1.0.x/LTS does not accept the 'package' argument.
    """
    from astropy.utils.data import get_pkg_data_filename
    try:
        f = get_pkg_data_filename(filename, package=package)
    except TypeError:
        from pkg_resources import resource_filename
        f = resource_filename(package, filename)
    return f
