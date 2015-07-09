#
# -*- coding: utf-8 -*-
from __future__ import print_function
#
#
def pca_qso(**kwargs):
    """Wrapper on pca_solve to handle QSO eigenspectra.
    """
    import os
    import os.path
    import pickle
    import pylab
    from astropy.io import ascii, fits
    import numpy as np
    from matplotlib.font_manager import fontManager, FontProperties
    from ...pydlutils.goddard.astro import get_juldate
    from ...pydlutils.image import djs_maskinterp
    from ...pydlutils.math import djs_median, computechi2
    from . import pca_solve, plot_eig, readspec, skymask, wavevector
    if 'inputfile' in kwargs:
        inputfile = kwargs['inputfile']
    else:
        inputfile = os.path.join(os.getenv('IDLSPEC2D_DIR'),
            'templates','eigeninput_qso.dat')
    if 'wavemin' in kwargs:
        wavemin = kwargs['wavemin']
    else:
        wavemin = 460.0
    if 'wavemax' in kwargs:
        wavemax = kwargs['wavemax']
    else:
        wavemax = 10000.0
    snmax = 100.0
    if 'niter' in kwargs:
        niter = kwargs['niter']
    else:
        niter = 200
    nkeep = 4
    minuse = 3
    #
    # Name the output files.
    #
    jd = get_juldate()
    outfile = "spEigenQSO-{0:d}".format(int(jd - 2400000.5))
    #
    # Read the input spectra
    #
    converters = {'plate': [ascii.convert_numpy(np.int32)],
        'mjd': [ascii.convert_numpy(np.int32)],
        'fiber': [ascii.convert_numpy(np.int32)] }
    input_data = ascii.read(inputfile,names=['plate','mjd','fiber','zfit'],converters=converters)
    plate = input_data['plate'].data
    mjd = input_data['mjd'].data
    fiber = input_data['fiber'].data
    zfit = input_data['zfit'].data
    spplate = readspec(plate,fiber,mjd=mjd,**kwargs)
    #
    # Insist that all of the requested spectra exist.
    #
    missing = spplate['plugmap']['FIBERID'] == 0
    if missing.any():
        imissing = missing.nonzero()[0]
        for k in imissing:
            print("Missing plate={0:d} mjd={1:d} fiber={2:d}".format(plate[k],mjd[k],fiber[k]))
        raise ValueError("{0:d} missing object(s).".format(missing.sum()))
    #
    # Do not fit where the spectrum may be dominated by sky-sub residuals.
    #
    objinvvar = skymask(spplate['invvar'],spplate['andmask'],spplate['ormask'])
    ifix = spplate['flux']**2 * objinvvar > snmax**2
    if ifix.any():
        objinvvar[ifix.nonzero()] = (snmax/spplate['flux'][ifix.nonzero()])**2
    #
    # Set the new wavelength mapping here.  If the binsz keyword is not set,
    # then bin size is determined from the first spectrum returned by readspec.
    # This is fine in the case where all spectra have the same bin size
    # (though their starting wavelengths may differ).  However, this may not
    # be a safe assumption in the future.
    #
    if 'binsz' in kwargs:
        objdloglam = kwargs['binsz']
    else:
        objdloglam = spplate['loglam'][0,1] - spplate['loglam'][0,0]
    newloglam = wavevector(np.log10(wavemin),np.log10(wavemax),
        binsz=objdloglam)
    #
    # Do PCA solution.
    #
    if 'qso' in kwargs:
        #
        # Solve for one component at a time, like the old pca_qso.pro
        #
        objflux = spplate['flux'].copy()
        objloglam = spplate['loglam'].copy()
        nobj, npix = spplate['flux'].shape
        for ikeep in range(nkeep):
            print("Solving for eigencomponent #{0:d} of {1:d}".format(ikeep+1,nkeep))
            pcaflux1 = pca_solve(objflux,objinvvar,objloglam,zfit,
                niter=niter,nkeep=1,newloglam=newloglam,aesthetics='mean')
            if ikeep == 0:
                #
                # Create new pcaflux dict
                #
                saveflux = pcaflux1['newflux']
                pcaflux = dict()
                for k in pcaflux1:
                    pcaflux[k] = pcaflux1[k].copy()
            else:
                #
                # Add to existing dict
                #
                # for k in pcaflux1:
                #     pcaflux[k] = np.vstack((pcaflux[k],pcaflux1[k]))
                pcaflux['flux'] = np.vstack((pcaflux['flux'],pcaflux1['flux']))
                pcaflux['eigenval'] = np.concatenate((pcaflux['eigenval'],pcaflux1['eigenval']))
            #
            # Re-solve for the coefficients using all PCA components so far
            #
            acoeff = np.zeros((nobj,ikeep+1),dtype=pcaflux1['acoeff'].dtype)
            for iobj in range(nobj):
                out = computechi2(saveflux[iobj,:],np.sqrt(pcaflux1['newivar'][iobj,:]),
                    pcaflux['flux'].T)
                acoeff[iobj,:] = out['acoeff']
            #
            # Prevent re-binning of spectra on subsequent calls to pca_solve()
            #
            objloglam = None
            if ikeep == 0:
                objflux = saveflux - np.outer(acoeff,pcaflux['flux'])
            else:
                objflux = saveflux - np.dot(acoeff,pcaflux['flux'])
            # objflux = saveflux - np.outer(acoeff,pcaflux['flux'])
            objinvvar = pcaflux1['newivar']
            pcaflux['acoeff'] = acoeff
    else:
        #
        # Do a normal simultaneous PCA solution
        #
        pcaflux = pca_solve(spplate['flux'],objinvvar,spplate['loglam'],zfit,
            niter=niter,nkeep=nkeep,newloglam=newloglam,aesthetics='mean')
    #
    # Fill in bad data with a running median of the good data.
    #
    qgood = pcaflux['usemask'] >= minuse
    medflux = np.zeros(pcaflux['flux'].shape,dtype=pcaflux['flux'].dtype)
    if not qgood.all():
        for i in range(nkeep):
            medflux[i,qgood] = djs_median(pcaflux['flux'][i,qgood],
                width=51,boundary='nearest')
            medflux[i,:] = djs_maskinterp(medflux[i,:],~qgood,const=True)
        pcaflux['flux'][:,~qgood] = medflux[:,~qgood]
    #
    # Dump input fluxes to a file for debugging purposes.
    #
    if 'dump' in kwargs:
        f = open(kwargs['dump'],'w')
        pickle.dump(pcaflux, f)
        f.close()
    #
    # Make plots
    #
    colorvec = ['k','r','g','b','m','c']
    smallfont = FontProperties(size='xx-small');
    nspectra = pcaflux['newflux'].shape[0]
    if 'flux' in kwargs:
        nfluxes = 30
        separation = 5.0
        nplots = nspectra/nfluxes
        if nspectra % nfluxes > 0:
            nplots += 1
        for k in range(nplots):
            istart = k*nfluxes
            iend = min(istart+nfluxes,nspectra) - 1
            fig = pylab.figure(dpi=100)
            ax = fig.add_subplot(111)
            for l in range(istart,iend+1):
                p = ax.plot(10.0**pcaflux['newloglam'],pcaflux['newflux'][l,:]+separation*(l%nfluxes),
                    '%s-'%colorvec[l%len(colorvec)],linewidth=1)
            ax.set_xlabel(r'Wavelength [$\AA$]')
            ax.set_ylabel(r'Flux [$\mathsf{10^{-17} erg\, cm^{-2} s^{-1} \AA^{-1}}$] + Constant')
            ax.set_title('QSOs: Input Spectra %4d-%4d' % (istart+1,iend+1))
            ax.set_ylim(pcaflux['newflux'][istart,:].min(),pcaflux['newflux'][iend-1,:].max()+separation*(nfluxes-1))
            fig.savefig('%s.flux.%04d-%04d.png'%(outfile,istart+1,iend+1))
            pylab.close(fig)
    fig = pylab.figure(dpi=100)
    ax = fig.add_subplot(111)
    p = ax.plot(10.0**pcaflux['newloglam'],(pcaflux['newivar'] == 0).sum(0)/float(nspectra),'k-')
    ax.set_xlabel(r'Wavelength [$\AA$]')
    ax.set_ylabel('Fraction of spectra with missing data')
    ax.set_title('Missing Data')
    fig.savefig(outfile+'.missing.png')
    pylab.close(fig)
    aratio10 = pcaflux['acoeff'][:,1]/pcaflux['acoeff'][:,0]
    aratio20 = pcaflux['acoeff'][:,2]/pcaflux['acoeff'][:,0]
    aratio30 = pcaflux['acoeff'][:,3]/pcaflux['acoeff'][:,0]
    fig = pylab.figure(dpi=100)
    ax = fig.add_subplot(111)
    p = ax.plot(aratio10,aratio20,marker='None',linestyle='None')
    for k in range(len(aratio10)):
        t = ax.text(aratio10[k],aratio20[k],'%04d-%04d'%(plate[k],fiber[k]),
            horizontalalignment='center', verticalalignment='center',
            color=colorvec[k%len(colorvec)],
            fontproperties=smallfont)
    # ax.set_xlim([aratio10.min(),aratio10.max])
    # ax.set_xlim([aratio20.min(),aratio20.max])
    ax.set_xlabel('Eigenvalue Ratio, $a_1/a_0$')
    ax.set_ylabel('Eigenvalue Ratio, $a_2/a_0$')
    ax.set_title('QSOs: Eigenvalue Ratios')
    fig.savefig(outfile+'.a2_v_a1.png')
    pylab.close(fig)
    fig = pylab.figure(dpi=100)
    ax = fig.add_subplot(111)
    p = ax.plot(aratio20,aratio30,marker='None',linestyle='None')
    for k in range(len(aratio10)):
        t = ax.text(aratio20[k],aratio30[k],'%04d-%04d'%(plate[k],fiber[k]),
            horizontalalignment='center', verticalalignment='center',
            color=colorvec[k%len(colorvec)],
            fontproperties=smallfont)
    # ax.set_xlim([aratio10.min(),aratio10.max])
    # ax.set_xlim([aratio20.min(),aratio20.max])
    ax.set_xlabel('Eigenvalue Ratio, $a_2/a_0$')
    ax.set_ylabel('Eigenvalue Ratio, $a_3/a_0$')
    ax.set_title('QSOs: Eigenvalue Ratios')
    fig.savefig(outfile+'.a3_v_a2.png')
    pylab.close(fig)
    #
    # Save output to FITS file.
    #
    if os.path.exists(outfile+'.fits'):
        os.remove(outfile+'.fits')
    hdu0 = fits.PrimaryHDU(pcaflux['flux'])
    hdu1 = fits.new_table(fits.ColDefs([
        fits.Column(name='plate',format='J',array=plate),
        fits.Column(name='mjd',format='J',array=mjd),
        fits.Column(name='fiber',format='J',array=fiber),
        fits.Column(name='redshift',format='D',array=zfit)]))
    hdulist = fits.HDUList([hdu0,hdu1])
    hdulist[0].header.update('OBJECT','QSO')
    hdulist[0].header.update('COEFF0',pcaflux['newloglam'][0])
    hdulist[0].header.update('COEFF1',pcaflux['newloglam'][1]-pcaflux['newloglam'][0])
    hdulist[0].header.update('IDLUTILS','pydlutils','Version of idlutils')
    hdulist[0].header.update('SPEC2D','eigenspectra','Version of idlspec2d')
    hdulist[0].header.update('RUN2D',os.getenv('RUN2D'),'Version of 2d reduction')
    hdulist[0].header.update('RUN1D',os.getenv('RUN1D'),'Version of 1d reduction')
    for i in range(len(pcaflux['eigenval'])):
        hdulist[0].header.update("EIGEN%d" % i,pcaflux['eigenval'][i])
    hdulist[1].header.update('FILENAME',inputfile)
    hdulist.writeto(outfile+'.fits')
    plot_eig(outfile+'.fits')
    return
