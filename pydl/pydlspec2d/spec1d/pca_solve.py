# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-


def pca_solve(flux, ivar, loglam=None, zfit=None, aesthetics='mean',
              newloglam=None, wavemin=None, wavemax=None,
              maxiter=0, niter=10, nkeep=3, nreturn=None, verbose=False):
    """Replacement for idlspec2d pca_solve.pro.

    Parameters
    ----------
    flux : array-like
        The input spectral flux.
    ivar : array-like
        The inverse variance of the spectral flux.
    loglam : array-like, optional
        The input wavelength solution.
    zfit : array-like, optional
        The redshift of each input spectrum.
    aesthetics : :class:`str`, optional
        This parameter will be passed to
        :func:`~pydl.pydlspec2d.spec2d.combine1fiber`.
    newloglam : array-like, optional
        The output wavelength solution.
    wavemin : :class:`float`, optional
        Minimum wavelength if `newloglam` is not specified.
    wavemax : :class:`float`, optional
        Maximum wavelength if `newloglam` is not specified.
    maxiter : :class:`int`, optional
        Stop PCA+reject iterations after this number.
    niter : :class:`int`, optional
        Stop PCA iterations after this number.
    nkeep : :class:`int`, optional
        Number of PCA components to keep.
    nreturn : :class:`int`, optional
        Number of PCA components to return, usually the same as `nkeep`.
    verbose : :class:`bool`, optional
        If ``True``, print extra information.

    Returns
    -------
    pca_solve : :class:`dict`
        The PCA solution.
    """
    import time
    import numpy as np
    from astropy import log
    from . import wavevector
    from ..spec2d import combine1fiber
    from ... import pcomp
    from ...pydlutils.math import computechi2, djs_reject
    if verbose:
        log.setLevel('DEBUG')
    if nreturn is None:
        nreturn = nkeep
    if len(flux.shape) == 1:
        nobj = 1
        npix = flux.shape[0]
    else:
        nobj, npix = flux.shape
    log.info("Building PCA from {0:d} object spectra.".format(nobj))
    #
    # The redshift of each object in pixels would be logshift/objdloglam.
    #
    if zfit is None:
        logshift = np.zeros((nobj,), dtype=flux.dtype)
    else:
        logshift = np.log10(1.0 + zfit)
    #
    # Determine the new wavelength mapping.
    #
    if loglam is None:
        newflux = flux
        newivar = ivar
        if newloglam is None:
            raise ValueError("newloglam must be set if loglam is not!")
        nnew = flux.shape[1]
    else:
        if newloglam is None:
            igood = loglam != 0
            dloglam = loglam[1] - loglam[0]
            logmin = loglam[igood].min() - logshift.max()
            logmax = loglam[igood].max() - logshift.min()
            if wavemin is not None:
                logmin = max(logmin, np.log10(wavemin))
            if wavemax is not None:
                logmax = min(logmax, np.log10(wavemax))
            fullloglam = wavevector(logmin, logmax, binsz=dloglam)
        else:
            fullloglam = newloglam
            dloglam = fullloglam[1] - fullloglam[0]
        nnew = fullloglam.size
        fullflux = np.zeros((nobj, nnew), dtype='d')
        fullivar = np.zeros((nobj, nnew), dtype='d')
        #
        # Shift each spectrum to z = 0 and sample at the output wavelengths
        #
        if loglam.ndim == 1:
            indx = loglam > 0
            rowloglam = loglam[indx]
        for iobj in range(nobj):
            log.info("OBJECT {0:5d}".format(iobj))
            if loglam.ndim > 1:
                if loglam.shape[0] != nobj:
                    raise ValueError('Wrong number of dimensions for loglam.')
                indx = loglam[iobj, :] > 0
                rowloglam = loglam[iobj, indx]
            flux1, ivar1 = combine1fiber(rowloglam-logshift[iobj],
                flux[iobj, indx], fullloglam, objivar=ivar[iobj, indx],
                binsz=dloglam, aesthetics=aesthetics, verbose=verbose)
            fullflux[iobj, :] = flux1
            fullivar[iobj, :] = ivar1
        #
        # Find the columns out side of which there is no data at all
        #
        # nzi = fullivar.nonzero()
        # firstcol = nzi[1].min()
        # lastcol = nzi[1].max()
        # newflux = fullflux[:, firstcol:lastcol+1]
        # newivar = fullivar[:, firstcol:lastcol+1]
        # newloglam = fullloglam[firstcol:lastcol+1]
        # nnew = newloglam.size
        newflux = fullflux
        newivar = fullivar
        newloglam = fullloglam
    nzi = newivar.nonzero()
    first_nonzero = (np.arange(nobj, dtype=nzi[0].dtype),
        np.array([nzi[1][nzi[0] == k].min() for k in range(nobj)]))
    #
    # Construct the synthetic weight vector, to be used when replacing the
    # low-S/N object pixels with the reconstructions.
    #
    synwvec = np.ones((nnew,), dtype='d')
    for ipix in range(nnew):
        indx = newivar[:, ipix] != 0
        if indx.any():
            synwvec[ipix] = newivar[indx, ipix].mean()
    fluxdict = {'flux': newflux, 'eigenval': 1.0, 'acoeff': 1.0,
        'outmask': np.ones((nnew,), dtype='i4'),
        'usemask': np.ones((nnew,), dtype='i4'),
        'newflux': newflux, 'newivar': newivar,
        'newloglam': newloglam, }   # 'emevecs': 1.0, 'emevals': 1.0}
    #
    # If there is only one object spectrum, then all we can do is return it.
    #
    if nobj == 1:
        fluxdict['flux'] = newflux.astype('f')
        return fluxdict
    #
    # Rejection iteration loop.
    #
    qdone = 0
    iiter = 0
    #
    # Begin with all points good.
    #
    outmask = None
    inmask = newivar != 0
    ymodel = None
    # emevecs, emevals = pydlutils.empca(newflux, inmask)
    # fluxdict['emevecs'] = emevecs
    # fluxdict['emevals'] = emeveals
    while qdone == 0 and iiter <= maxiter:
        log.debug('starting djs_reject')
        outmask, qdone = djs_reject(newflux, ymodel, inmask=inmask,
                                    outmask=outmask, invvar=newivar)
        log.debug('finished with djs_reject')
        #
        # Iteratively do the PCA solution
        #
        filtflux = newflux.copy()
        acoeff = np.zeros((nobj, nkeep), dtype='d')
        t0 = time.time()
        for ipiter in range(niter):
            #
            # We want to get these values from the pcomp routine.
            #
            # eigenval = 1
            # coeff = 1
            flux0 = np.tile(filtflux[first_nonzero], nnew).reshape(nnew, nobj).transpose()
            # flux0 = np.tile(filtflux, nnew).reshape(nnew, nobj).transpose()
            totflux = np.absolute(filtflux - flux0).sum(1)
            goodobj = totflux > 0
            if goodobj.all():
                tmp = pcomp(filtflux.T)  # , standardize=True)
                pres = tmp.derived
                eigenval = tmp.eigenvalues
            else:
                tmp = pcomp(filtflux[goodobj, :].T)  # , standardize=True)
                pres = np.zeros((nobj, nnew), dtype='d')
                pres[goodobj, :] = tmp.derived
                eigenval = np.zeros((nobj,), dtype='d')
                eigenval[goodobj] = tmp.eigenvalues
            maskivar = newivar * outmask
            sqivar = np.sqrt(maskivar)
            for iobj in range(nobj):
                out = computechi2(newflux[iobj, :], sqivar[iobj, :],
                    pres[:, 0:nkeep])
                filtflux[iobj, :] = (maskivar[iobj, :] * newflux[iobj, :] +
                    synwvec*out.yfit) / (maskivar[iobj, :] + synwvec)
                acoeff[iobj, :] = out.acoeff
            log.info("The elapsed time for iteration #{0:2d} is {1:6.2f} s.".format(ipiter+1, time.time()-t0))
        #
        # Now set ymodel for rejecting points.
        #
        ymodel = np.dot(acoeff, pres[:, 0:nkeep].T)
        iiter += 1
    if nobj == 1:
        usemask = outmask
    else:
        usemask = outmask.sum(0)
    fluxdict['usemask'] = usemask
    fluxdict['outmask'] = outmask
    fluxdict['flux'] = pres[:, 0:nreturn].transpose().astype('f')
    fluxdict['eigenval'] = eigenval[0:nreturn]
    fluxdict['acoeff'] = acoeff
    return fluxdict
