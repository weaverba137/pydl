# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-


def pca_solve(newflux, newivar, newloglam,
              maxiter=0, niter=10, nkeep=3, nreturn=None, verbose=False):
    """Replacement for idlspec2d pca_solve.pro.

    Parameters
    ----------
    newflux : array-like
        The input spectral flux, assumed to have a common wavelength and
        redshift system.
    newivar : array-like
        The inverse variance of the spectral flux.
    newloglam : array-like, optional
        The output wavelength solution.
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
    :class:`dict`
        The PCA solution.
    """
    import time
    import numpy as np
    from astropy import log
    from ... import pcomp
    from ...pydlutils.math import computechi2, djs_reject
    if verbose:
        log.setLevel('DEBUG')
    if nreturn is None:
        nreturn = nkeep
    if len(newflux.shape) == 1:
        nobj = 1
        npix = newflux.shape[0]
    else:
        nobj, npix = newflux.shape
    log.info("Building PCA from {0:d} object spectra.".format(nobj))
    nnew = newloglam.size
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
    fluxdict = dict()
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
