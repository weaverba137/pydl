# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from __future__ import print_function
#
def pca_solve(flux,ivar,loglam=None,zfit=None,**kwargs):
    """Replacement for idlspec2d pca_solve.pro
    """
    import time
    import numpy as np
    import pydl.pydlspec2d.spec1d
    import pydl.pydlspec2d.spec2d
    from ... import pcomp
    from ...pydlutils.math import computechi2, djs_reject
    if 'maxiter' in kwargs:
        maxiter = kwargs['maxiter']
    else:
        maxiter = 0
    if 'niter' in kwargs:
        niter = kwargs['niter']
    else:
        niter = 10
    if 'nkeep' in kwargs:
        nkeep = kwargs['nkeep']
    else:
        nkeep = 3
    if 'nreturn' in kwargs:
        nreturn = kwargs['nreturn']
    else:
        nreturn = nkeep
    if len(flux.shape) == 1:
        nobj = 1
        npix = flux.shape[0]
    else:
        nobj, npix = flux.shape
    print("Building PCA from {0:d} object spectra.".format(nobj))
    #
    # The redshift of each object in pixels would be logshift/objdloglam.
    #
    if zfit is None:
        logshift = np.zeros((nobj,),dtype=flux.dtype)
    else:
        logshift = np.log10(1.0 + zfit)
    #
    # Determine the new wavelength mapping.
    #
    if loglam is None:
        newflux = flux
        newivar = ivar
        newloglam = kwargs['newloglam']
        nnew = flux.shape[1]
    else:
        if 'newloglam' in kwargs:
            fullloglam = kwargs['newloglam']
            dloglam = fullloglam[1] - fullloglam[0]
        else:
            igood = loglam != 0
            dloglam = loglam[1] - loglam[0]
            logmin = loglam[igood].min() - logshift.max()
            logmax = loglam[igood].max() - logshift.min()
            if 'wavemin' in kwargs:
                logmin = max(logmin, np.log10(wavemin))
            if 'wavemax' in kwargs:
                logmax = min(logmax,np.log10(wavemax))
            fullloglam = pydl.pydlspec2d.spec1d.wavevector(logmin,logmax,binsz=dloglam)
        nnew = fullloglam.size
        fullflux = np.zeros((nobj,nnew),dtype='d')
        fullivar = np.zeros((nobj,nnew),dtype='d')
        #
        # Shift each spectrum to z = 0 and sample at the output wavelengths
        #
        if loglam.ndim == 1:
            indx = loglam > 0
            rowloglam = loglam[indx]
        for iobj in range(nobj):
            print("OBJECT {0:5d}".format(iobj))
            if loglam.ndim > 1:
                if loglam.shape[0] != nobj:
                    raise ValueError('Wrong number of dimensions for loglam.')
                indx = loglam[iobj,:] > 0
                rowloglam = loglam[iobj,indx]
            flux1,ivar1 = pydl.pydlspec2d.spec2d.combine1fiber(rowloglam-logshift[iobj],flux[iobj,indx],
                ivar[iobj,indx],newloglam=fullloglam,binsz=dloglam,aesthetics='mean') # ,verbose=True)
            fullflux[iobj,:] = flux1
            fullivar[iobj,:] = ivar1
        #
        # Find the columns out side of which there is no data at all
        #
        # nzi = fullivar.nonzero()
        # firstcol = nzi[1].min()
        # lastcol = nzi[1].max()
        # newflux = fullflux[:,firstcol:lastcol+1]
        # newivar = fullivar[:,firstcol:lastcol+1]
        # newloglam = fullloglam[firstcol:lastcol+1]
        # nnew = newloglam.size
        newflux = fullflux
        newivar = fullivar
        newloglam = fullloglam
    nzi = newivar.nonzero()
    first_nonzero = (np.arange(nobj,dtype=nzi[0].dtype),
        np.array([nzi[1][nzi[0]==k].min() for k in range(nobj)]))
    #
    # Construct the synthetic weight vector, to be used when replacing the
    # low-S/N object pixels with the reconstructions.
    #
    synwvec = np.ones((nnew,),dtype='d')
    for ipix in range(nnew):
        indx = newivar[:,ipix] != 0
        if indx.any():
            synwvec[ipix] = newivar[indx,ipix].mean()
    fluxdict = {'flux':newflux, 'eigenval':1.0, 'acoeff':1.0,
        'outmask':np.ones((nnew,),dtype='i4'),
        'usemask':np.ones((nnew,),dtype='i4'),
        'newflux':newflux,'newivar':newivar,
        'newloglam':newloglam, }
        # 'emevecs':1.0, 'emevals':1.0}
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
    # emevecs, emevals = pydlutils.empca(newflux,inmask)
    # fluxdict['emevecs'] = emevecs
    # fluxdict['emevals'] = emeveals
    while qdone == 0 and iiter <= maxiter:
        # print('starting djs_reject')
        outmask,qdone = djs_reject(newflux,ymodel,
            inmask=inmask,outmask=outmask,
            invvar=newivar)
        # print('finished with djs_reject')
        #
        # Iteratively do the PCA solution
        #
        filtflux = newflux.copy()
        acoeff = np.zeros((nobj,nkeep),dtype='d')
        t0 = time.time()
        for ipiter in range(niter):
            #
            # We want to get these values from the pcomp routine.
            #
            # eigenval = 1
            # coeff = 1
            flux0 = np.tile(filtflux[first_nonzero],nnew).reshape(nnew,nobj).transpose()
            # flux0 = np.tile(filtflux,nnew).reshape(nnew,nobj).transpose()
            totflux = np.absolute(filtflux - flux0).sum(1)
            goodobj = totflux > 0
            if goodobj.all():
                tmp = pcomp(filtflux.T) # ,standardize=True)
                pres = tmp.derived
                eigenval = tmp.eigenvalues
            else:
                tmp = pcomp(filtflux[goodobj,:].T) # ,standardize=True)
                pres = np.zeros((nobj,nnew),dtype='d')
                pres[goodobj,:] = tmp.derived
                eigenval = np.zeros((nobj,),dtype='d')
                eigenval[goodobj] = tmp.eigenvalues
            maskivar = newivar * outmask
            sqivar = np.sqrt(maskivar)
            for iobj in range(nobj):
                out = computechi2(newflux[iobj,:], sqivar[iobj,:],
                    pres[:,0:nkeep])
                filtflux[iobj,:] = (maskivar[iobj,:] * newflux[iobj,:] +
                    synwvec*out.yfit) / (maskivar[iobj,:] + synwvec)
                acoeff[iobj,:] = out.acoeff
            if 'quiet' not in kwargs:
                print("The elapsed time for iteration #{0:2d} is {1:6.2f} s.".format(ipiter+1,time.time()-t0))
        #
        # Now set ymodel for rejecting points.
        #
        ymodel = np.dot(acoeff,pres[:,0:nkeep].T)
        iiter += 1
    if nobj == 1:
        usemask = outmask
    else:
        usemask = outmask.sum(0)
    fluxdict['usemask'] = usemask
    fluxdict['outmask'] = outmask
    fluxdict['flux'] = pres[:,0:nreturn].transpose().astype('f')
    fluxdict['eigenval'] = eigenval[0:nreturn]
    fluxdict['acoeff'] = acoeff
    return fluxdict

