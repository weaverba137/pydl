# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from __future__ import print_function
#
def combine1fiber(inloglam,objflux,objivar=None,**kwargs):
    """Combine several spectra of the same object, or resample a single spectrum.
    """
    import numpy as np
    import pydl.pydlutils.bspline
    import pydl.pydlutils.sdss
    from pydl import smooth
    from . import aesthetics
    #
    # Check that dimensions of inputs are valid.
    #
    if 'newloglam' in kwargs:
        newloglam = kwargs['newloglam']
        nfinalpix = len(newloglam)
    else:
        raise ValueError('newloglam is required.')
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
        npix = inloglam.shape[0]
        nspec = 1
        specnum = np.zeros(inloglam.shape,dtype=inloglam.dtype)
    else:
        nspec,npix = inloglam.shape
        specnum = np.tile(np.arange(nspec),npix).reshape(npix,nspec).transpose()
    #
    # Use fullcombmask for modifying the pixel masks in the original input files.
    #
    fullcombmask = np.zeros(npix)
    newflux = np.zeros(nfinalpix,dtype=inloglam.dtype)
    newmask = np.zeros(nfinalpix,dtype='i4')
    newivar = np.zeros(nfinalpix,dtype=inloglam.dtype)
    newdisp = np.zeros(nfinalpix,dtype=inloglam.dtype)
    newsky = np.zeros(nfinalpix,dtype=inloglam.dtype)
    newdispweight = np.zeros(nfinalpix,dtype=inloglam.dtype)
    if objivar is None:
        nonzero = np.arange(npix,dtype='i4')
        ngood = npix
    else:
        nonzero = (objivar > 0).nonzero()[0]
        ngood = len(nonzero)
    #
    # ormask is needed to create andmask
    #
    andmask = np.zeros(nfinalpix,dtype='i4')
    ormask = np.zeros(nfinalpix,dtype='i4')
    if ngood == 0:
        #
        # In this case of no good points, set the nodata bit everywhere.
        # Also if noplug is set in the first input bit-mask, assume it
        # should be set everywhere in the output bit masks.  No other bits
        # are set.
        #
        if 'verbose' in kwargs:
            print('No good points')
        bitval = pydl.pydlutils.sdss.sdss_flagval('SPPIXMASK','NODATA')
        if 'finalmask' in kwargs:
            bitval = bitval | (pydl.pydlutils.sdss.sdss_flagval('SPPIXMASK','NOPLUG') *
                (finalmask[0] & pydl.pydlutils.sdss.sdss_flagval('SPPIXMASK','NODATA')))
        andmask = andmask | bitval
        ormask = ormask | bitval
        return (newflux,newivar)
    else:
        #
        # Now let's break sorted wavelengths into groups where pixel
        # separations are larger than maxsep.
        #
        isort = nonzero[inloglam[nonzero].argsort()]
        wavesort = inloglam[isort]
        padwave = np.insert(wavesort,0,wavesort.min() - 2.0*maxsep)
        padwave = np.append(padwave,wavesort.max() + 2.0*maxsep)
        ig1 = ((padwave[1:ngood+1]-padwave[0:ngood]) > maxsep).nonzero()[0]
        ig2 = ((padwave[2:ngood+2]-padwave[1:ngood+1]) > maxsep).nonzero()[0]
        if ig1.size != ig2.size:
            raise ValueError('Grouping tricks did not work!')
        for igrp in range(ig1.size):
            ss = isort[ig1[igrp]:ig2[igrp]+1]
            if ss.size > 2:
                if objivar is None:
                    #
                    # Fit without variance
                    #
                    sset,bmask = pydl.pydlutils.bspline.iterfit(inloglam[ss],objflux[ss],
                        nord=nord,groupbadpix=True,requiren=1,bkspace=bkptbin,
                        silent=True)
                else:
                    #
                    # Fit with variance
                    #
                    sset,bmask = pydl.pydlutils.bspline.iterfit(inloglam[ss],objflux[ss],
                        invvar=objivar[ss],
                        nord=nord,groupbadpix=True,requiren=1,bkspace=bkptbin,
                        silent=True)
                if np.sum(np.absolute(sset.coeff)) == 0:
                    sset = None
                    bmask = np.zeros(len(ss))
                    if 'verbose' in kwargs:
                        print('WARNING: All B-spline coefficients have been set to zero!')
            else:
                bmask = np.zeros(len(ss))
                sset = None
                if 'verbose' in kwargs:
                    print('WARNING: Not enough data for B-spline fit!')
            inside = ((newloglam >= inloglam[ss].min()-EPS) &
                (newloglam <= inloglam[ss].max()+EPS)).nonzero()[0]
            #
            # It is possible for numinside to be zero, if the input data points
            # span an extremely small wavelength range, within which there are
            # no output wavelengths.
            #
            if sset is not None and len(inside) > 0:
                newflux[inside],bvalumask = sset.value(newloglam[inside])
                if bvalumask.any():
                    newmask[inside[bvalumask]] = 1
                if 'verbose' in kwargs:
                    print('Masked {0:d} of {1:d} pixels.'.format(bmask.sum()-bmask.size,bmask.size))
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
                        objivar[ss[ireplace]] = 0.0
                        if 'verbose' in kwargs:
                            print('Replaced {0:d} pixels in objivar.'.format(len(ss[ireplace])))
                    if 'finalmask' in kwargs:
                        finalmask[ss[ireplace]] = (finalmask[ss[ireplace]] |
                            pydl.pydlutils.sdss.sdss_flagval('SPPIXMASK','COMBINEREJ'))
            fullcombmask[ss] = bmask
        #
        # Combine inverse variance and pixel masks.
        #
        # Start with all bits set in andmask
        #
        andmask[:] = -1
        for j in range(int(specnum.max())+1):
            these = specnum == j
            if these.any():
                inbetween = ((newloglam >= inloglam[these].min()) &
                    (newloglam <= inloglam[these].max()))
                if inbetween.any():
                    jnbetween = inbetween.nonzero()[0]
                    #
                    # Conserve inverse variance by doing a linear interpolation
                    # on that quantity.
                    #
                    result = np.interp(newloglam[jnbetween],
                        inloglam[these],objivar[these]*fullcombmask[these])
                    #
                    # Grow the fullcombmask below to reject any new sampling
                    # containing even a partial masked pixel.
                    #
                    smask = np.interp(newloglam[jnbetween],inloglam[these],
                        fullcombmask[these].astype(inloglam.dtype))
                    result *= smask >= (1.0 - EPS)
                    newivar[jnbetween] += result*newmask[jnbetween]
                lowside = np.floor((inloglam[these]-newloglam[0])/binsz).astype('i4')
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
                    newdisp[jnbetween] += result * np.interp(newloglam[jnbetween],
                        inloglam[these],indisp[these])
                    newsky[jnbetween] += result * np.interp(newloglam[jnbetween],
                        inloglam[these],skyflux[these])
        if 'indisp' in kwargs:
            newdisp /= newdispweight + (newdispweight == 0)
            newsky /= newdispweight + (newdispweight == 0)
    #
    # Grow regions where 3 or more pixels are rejected together ???
    #
    # print(newivar)
    foo = smooth(newivar,3)
    # print(foo)
    # sys.exit(1)
    badregion = np.absolute(foo) < EPS
    if badregion.any():
        if 'verbose' in kwargs:
            print('WARNING: Growing bad pixel region, {0:d} pixels found.'.format(badregion.sum()))
        ibad = badregion.nonzero()[0]
        lowerregion = np.where(ibad-2 < 0,0,ibad-2)
        upperregion = np.where(ibad+2 > nfinalpix-1,nfinalpix-1,ibad+2)
        newivar[lowerregion] = 0.0
        newivar[upperregion] = 0.0
    #
    # Replace NaNs in combined spectra; this should really never happen.
    #
    inff = ((~np.isfinite(newflux)) | (~np.isfinite(newivar)))
    if inff.any():
        print('WARNING: {0:d} NaNs in combined spectra.'.format(inff.sum()))
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
    newflux = aesthetics(newflux,newivar,method=amethod)
    # if 'interpolate' in kwargs:
    #     newflux = pydlutils.image.djs_maskinterp(newflux,~goodpts,const=True)
    # else:
    #     newflux[~goodpts] = newflux[goodpts].mean()
    if goodpts.any():
        minglam = newloglam[goodpts].min()
        maxglam = newloglam[goodpts].max()
        ibad = ((newloglam < minglam) | (newloglam > maxglam))
        if ibad.any():
            ormask[ibad] |= pydl.pydlutils.sdss.sdss_flagval('SPPIXMASK','NODATA')
    #
    # Replace values of -1 in the andmask with 0.
    #
    andmask *= (andmask != -1)
    #
    # Copy the nodata bad pixels from the ormask to the andmask.
    #
    andmask |= ormask & pydl.pydlutils.sdss.sdss_flagval('SPPIXMASK','NODATA')
    return (newflux,newivar)

