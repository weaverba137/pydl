# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from __future__ import print_function
#
def readspec(platein,fiber='all',**kwargs):
    """Read SDSS/BOSS spec2d & spec1d files.

    INPUTS:
        plate - Plate number(s)
    OPTIONAL INPUTS:
        topdir - Override the value of $BOSS_SPECTRO_REDUX.
        run2d  - Override the value of $RUN2D.
        run1d  - Override the value of $RUN1D
        path   - Override all path information with this directory name.
        align  - If set, align all the spectra in wavelength.
        znum   - If set, return the znum-th best fit reshift fit, instead of the best
    """
    import os
    import os.path
    from astropy.io import fits as pyfits
    import numpy as np
    from . import number_of_fibers, latest_mjd, spec_path, spec_append
    try:
        nplate = len(platein)
        plate = platein
    except TypeError:
        nplate = 1
        plate = np.array([platein],dtype='i4')
    if 'mjd' in kwargs:
        try:
            nmjd = len(kwargs['mjd'])
        except TypeError:
            nmjd = 1
        if nmjd != nplate:
            raise TypeError("Plate & MJD must have the same length!")
    if 'run2d' in kwargs:
        run2d = kwargs['run2d']
    else:
        run2d = os.getenv('RUN2D')
    if 'run1d' in kwargs:
        run1d = kwargs['run1d']
    else:
        run1d = os.getenv('RUN1D')
    if fiber == 'all':
        #
        # Read all fibers
        #
        nfibers = number_of_fibers(plate,**kwargs)
        total_fibers = nfibers.sum()
        platevec = np.zeros(total_fibers,dtype='i4')
        fibervec = np.zeros(total_fibers,dtype='i4')
        k = 0
        for p in np.unique(plate):
            n = np.unique(nfibers[plate == p])[0]
            platevec[k:k+n] = p
            fibervec[k:k+n] = np.arange(n) + 1
            k += n
    else:
        try:
            nfiber = len(fiber)
        except TypeError:
            nfiber = 1
        if nplate > 1 and nfiber > 1 and nplate != nfiber:
            raise TypeError("Plate & Fiber must have the same length!")
        if nplate > 1:
            platevec = np.array(plate,dtype='i4')
        else:
            platevec = np.zeros(nfiber,dtype='i4') + plate
        if nfiber > 1:
            fibervec = np.array(fiber,dtype='i4')
        else:
            fibervec = np.zeros(nplate,dtype='i4') + fiber
    if 'mjd' in kwargs:
        mjdvec = np.zeros(nplate,dtype='i4') + kwargs['mjd']
    else:
        mjdvec = latest_mjd(platevec,**kwargs)
    #
    # Now select unique plate-mjd combinations & read them
    #
    pmjd = ((np.array(platevec,dtype='u8') << 16) +
        np.array(mjdvec,dtype='u8'))
    #print(pmjd)
    upmjd = np.unique(pmjd)
    zupmjd = list(zip(upmjd>>16,upmjd&((1<<16)-1)))
    #print(zupmjd)
    spplate_data = dict()
    hdunames = ('flux','invvar','andmask','ormask','disp','plugmap','sky','loglam',)
    for thisplate,thismjd in zupmjd:
        #thisplate = int(p>>16)
        #thismjd = int(np.bitwise_and(p,(1<<16)-1))
        pmjdindex = ((platevec==thisplate) & (mjdvec==thismjd)).nonzero()[0]
        thisfiber = fibervec[pmjdindex]
        #print(type(thisplate),type(thismjd))
        pmjdstr = "{0:04d}-{1:05d}".format(int(thisplate),int(thismjd))
        if 'path' in kwargs:
            sppath = [ kwargs['path'] ]
        else:
            sppath = spec_path(thisplate,run2d=run2d)
        spfile = os.path.join(sppath[0],"spPlate-{0}.fits".format(pmjdstr))
        print(spfile)
        spplate = pyfits.open(spfile)
        #
        # Get wavelength coefficients from primary header
        #
        npix = spplate[0].header['NAXIS1']
        c0 = spplate[0].header['COEFF0']
        c1 = spplate[0].header['COEFF1']
        coeff0 = np.zeros(thisfiber.size,dtype='d') + c0
        coeff1 = np.zeros(thisfiber.size,dtype='d') + c1
        loglam0 = c0 + c1*np.arange(npix,dtype='d')
        loglam = np.resize(loglam0,(thisfiber.size,npix))
        #
        # Read the data images
        #
        for k in range(len(hdunames)):
            try:
                tmp = spplate[k].data[thisfiber-1,:]
            except IndexError:
                tmp = loglam
            if hdunames[k] not in spplate_data:
                if k == 0:
                    allpmjdindex = pmjdindex
                    allcoeff0 = coeff0
                    allcoeff1 = coeff1
                #
                # Put the data into the return structure
                #
                if hdunames[k] == 'plugmap':
                    spplate_data['plugmap'] = dict()
                    for c in spplate[5].columns.names:
                        spplate_data['plugmap'][c] = tmp.field(c)
                else:
                    spplate_data[hdunames[k]] = tmp
            else:
                #
                # Append data
                #
                if k == 0:
                    allpmjdindex = np.concatenate((allpmjdindex,pmjdindex))
                    if 'align' in kwargs:
                        mincoeff0 = min(allcoeff0)
                        if mincoeff0 == 0 and coeff0[0] > 0:
                            allcoeff0 = coeff0[0]
                            allcoeff1 = coeff1[1]
                        if mincoeff0 > 0 and coeff0[0] == 0:
                            coeff0 = mincoeff0
                            coeff1 = allcoeff1[0]
                        ps = np.floor( (coeff0[0] - mincoeff0)/coeff1[0] + 0.5)
                        if ps > 0:
                            coeff0 = coeff0 - ps*coeff1
                        else:
                            allcoeff0 = allcoeff0 + ps*allcoeff1
                    else:
                        ps = 0
                    allcoeff0 = np.concatenate((allcoeff0,coeff0))
                    allcoeff1 = np.concatenate((allcoeff1,coeff1))
                if hdunames[k] == 'plugmap':
                    for c in spplate[5].columns.names:
                        spplate_data['plugmap'][c] = np.concatenate(
                            (spplate_data['plugmap'][c],tmp.field(c)))
                else:
                    spplate_data[hdunames[k]] = spec_append(spplate_data[hdunames[k]],tmp,pixshift=ps)
        spplate.close()
        #
        # Read photoPlate information, if available
        #
        photofile = os.path.join(sppath[0],"photoPlate-{0}.fits".format(pmjdstr))
        if not os.path.exists(photofile):
            #
            # Hmm, maybe this is an SDSS-I,II plate
            #
            photofile = os.path.join(os.getenv('SPECTRO_MATCH'),run2d,
                os.path.basename(os.getenv('PHOTO_RESOLVE')),"{0:04d}".format(int(thisplate)),
                "photoPlate-{0}.fits".format(pmjdstr))
        if os.path.exists(photofile):
            photop = pyfits.open(photofile)
            tmp = photop[1].data[thisfiber-1]
            if 'tsobj' not in spplate_data:
                spplate_data['tsobj'] = dict()
                for c in photop[1].columns.names:
                    spplate_data['tsobj'][c] = tmp.field(c)
            else:
                for c in photop[1].columns.names:
                    spplate_data['tsobj'][c] = np.concatenate(
                        (spplate_data['tsobj'][c],tmp.field(c)))
            photop.close()

        #
        # Read redshift information, if available.
        #
        if 'znum' in kwargs:
            zfile = os.path.join(sppath[0],run1d,"spZall-{0}.fits".format(pmjdstr))
        else:
            zfile = os.path.join(sppath[0],run1d,"spZbest-{0}.fits".format(pmjdstr))
        if os.path.exists(zfile):
            spz = pyfits.open(zfile)
            if 'znum' in kwargs:
                nper = spz[0].header['DIMS0']
                zfiber = (thisfiber-1)*nper + kwargs['znum'] - 1
            else:
                zfiber = thisfiber
            tmp = spz[1].data[zfiber-1]
            if 'zans' not in spplate_data:
                spplate_data['zans'] = dict()
                for c in spz[1].columns.names:
                    spplate_data['zans'][c] = tmp.field(c)
            else:
                for c in spz[1].columns.names:
                    spplate_data['zans'][c] = np.concatenate(
                        (spplate_data['zans'][c],tmp.field(c)))
            spz.close()
    #
    # Reorder the data.  At this point allpmjdindex is an index for which
    # fiber[allpmjdindex] == spplate['plugmap']['FIBERID'], so we have to
    # reverse this mapping.
    #
    j = allpmjdindex.argsort()
    for k in spplate_data:
        if isinstance(spplate_data[k],dict):
            for c in spplate_data[k]:
                if spplate_data[k][c].ndim == 2:
                    spplate_data[k][c] = spplate_data[k][c][j,:]
                else:
                    spplate_data[k][c] = spplate_data[k][c][j]
        else:
            spplate_data[k] = spplate_data[k][j,:]
    allcoeff0 = allcoeff0[j]
    allcoeff1 = allcoeff1[j]
    #
    # If necessary, recompute the wavelengths
    #
    nfibers,npixmax = spplate_data['flux'].shape
    if 'align' in kwargs:
        loglam0 = allcoeff0[0] + allcoeff1[1]*np.arange(npixmax,dtype='d')
        spplate_data['loglam'] = np.resize(loglam0,(nfibers,npixmax))
    return spplate_data

