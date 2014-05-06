#
# -*- coding: utf-8 -*-
from __future__ import print_function
#
#
def pca_star(**kwargs):
    """Wrapper on pca_solve to handle stellar eigenspectra.
    """
    import os
    import os.path
    import pickle
    import pylab
    from astropy.io import fits as pyfits
    import numpy as np
    import pydl.pydlutils.yanny as yanny
    from matplotlib.font_manager import fontManager, FontProperties
    from ... import uniq
    from ...pydlutils.goddard.astro import get_juldate
    from ...pydlutils.image import djs_maskinterp
    from . import pca_solve, readspec, skymask
    if 'inputfile' in kwargs:
        inputfile = kwargs['inputfile']
    else:
        inputfile = os.path.join(os.getenv('IDLSPEC2D_DIR'),
            'templates','eigeninput_star.par')
    wavemin = 0
    wavemax = 0
    snmax = 100.0
    if 'niter' in kwargs:
        niter = kwargs['niter']
    else:
        niter = 10
    cspeed = 2.99792458e5
    #
    # Name the output files.
    #
    jd = get_juldate()
    outfile = "spEigenStar-{0:d}".format(int(jd - 2400000.5))
    #
    # Read the input spectra
    #
    par = yanny.yanny(inputfile,np=True)
    slist = par['EIGENOBJ']
    spplate = readspec(slist['plate'],slist['fiberid'],mjd=slist['mjd'],align=True,**kwargs)
    objdloglam = spplate['loglam'][0,1] - spplate['loglam'][0,0]
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
    # Find the list of unique star types
    #
    isort = np.argsort(slist['class'])
    classlist = slist['class'][isort[uniq(slist['class'][isort])]]
    #
    # Loop over each star type
    #
    fullflux = None
    namearr = list()
    colorvec = ['k','r','g','b','m','c']
    smallfont = FontProperties(size='xx-small');
    for c in classlist:
        #
        # Find the subclasses for this stellar type
        #
        print("Finding eigenspectra for Stellar class {0}".format(c))
        indx = (slist['class'] == c).nonzero()[0]
        nindx = indx.size
        thesesubclass = slist['subclass'][indx]
        isort = np.argsort(thesesubclass)
        subclasslist = thesesubclass[isort[uniq(thesesubclass[isort])]]
        nsubclass = subclasslist.size
        #
        # Solve for 2 eigencomponents if we have specified subclasses for
        # this stellar type
        #
        if nsubclass == 1:
            nkeep = 1
        else:
            nkeep = 2
        newloglam = spplate['loglam'][0,:]
        pcaflux = pca_solve(spplate['flux'][indx,:],objinvvar[indx,:],
            spplate['loglam'][indx,:], slist['cz'][indx]/cspeed, wavemin=wavemin,
            wavemax=wavemax, niter=niter, nkeep=nkeep, newloglam=newloglam)
        #
        # Interpolate over bad flux values in the middle of a spectrum,
        # and set fluxes to zero at the blue+red ends of the spectrum
        #
        # minuse = 1 # ?
        minuse = np.floor((nindx+1) / 3.0)
        qbad = pcaflux['usemask'] < minuse
        #
        # Interpolate over all bad pixels
        #
        for j in range(nkeep):
            pcaflux['flux'][j,:] = djs_maskinterp(pcaflux['flux'][j,:],qbad,const=True)
        #
        # Set bad pixels at the very start or end of the spectrum to zero instead
        #
        npix = qbad.size
        igood = (~qbad).nonzero()[0]
        if qbad[0]:
            pcaflux['flux'][:,0:igood[0]-1] = 0
        if qbad[npix-1]:
            pcaflux['flux'][:,igood[::-1][0]+1:npix] = 0
        #
        # Re-normalize the first eigenspectrum to a mean of 1
        #
        # print(pcaflux['flux'].dtype)
        norm = pcaflux['flux'][0,:].mean()
        pcaflux['flux'] /= norm
        pcaflux['acoeff'] *= norm
        # print(pcaflux['flux'].dtype)
        #
        # Now loop through each stellar subclass and reconstruct
        # an eigenspectrum for that subclass
        #
        thesesubclassnum = np.zeros(thesesubclass.size,dtype='i4')
        #
        # Create a figure for this class
        #
        fig = pylab.figure(dpi=100)
        ax = fig.add_subplot(111)
        for isub in range(nsubclass):
            ii = (thesesubclass == subclasslist[isub]).nonzero()[0]
            thesesubclassnum[ii] = isub
            if nkeep == 1:
                thisflux = pcaflux['flux'].reshape(pcaflux['newloglam'].shape)
                # print(thisflux.dtype)
            else:
                aratio = pcaflux['acoeff'][ii,1]/pcaflux['acoeff'][ii,0]
                #
                # np.median(foo) is equivalent to MEDIAN(foo,/EVEN)
                #
                thisratio = np.median(aratio)
                thisflux = pcaflux['flux'][0,:] + thisratio.astype('f') * pcaflux['flux'][1,:]
                # print(thisflux.dtype)
            if fullflux is None:
                fullflux = thisflux
            else:
                fullflux = np.vstack((fullflux,thisflux))
            namearr.append(subclasslist[isub])
            #
            # Plot spectra
            #
            plotflux = thisflux/thisflux.max()
            # print(pcaflux['newloglam'].shape)
            # print(plotflux.shape)
            ax.plot(10.0**pcaflux['newloglam'],plotflux,"{0}-".format(colorvec[isub%len(colorvec)]),linewidth=1)
            if isub == 0:
                ax.set_xlabel(r'Wavelength [$\AA$]')
                ax.set_ylabel('Flux [arbitrary units]')
                ax.set_title('STAR {0}: Eigenspectra Reconstructions'.format(c))
            t = ax.text(10.0**pcaflux['newloglam'][-1],plotflux[-1],subclasslist[isub],
                horizontalalignment='right',verticalalignment='center',
                color=colorvec[isub%len(colorvec)],fontproperties=smallfont)
        fig.savefig(outfile+'.{0}.png'.format(c))
        pylab.close(fig)
        #
        # Plot eigenvalue ratios if nkeep > 1
        #
        if nkeep > 1:
            fig = pylab.figure(dpi=100)
            ax = fig.add_subplot(111)
            allratio = pcaflux['acoeff'][:,1]/pcaflux['acoeff'][:,0]
            isort = thesesubclassnum.argsort()
            p = ax.plot(thesesubclassnum[isort],allratio[isort],marker='None',linestyle='None')
            for k in range(len(indx)):
                t = ax.text(thesesubclassnum[isort[k]],allratio[isort[k]],
                    "%04d-%04d"%(slist['plate'][indx[isort[k]]],slist['fiberid'][indx[isort[k]]]),
                    horizontalalignment='center', verticalalignment='center',
                    color=colorvec[k%len(colorvec)],
                    fontproperties=smallfont)
            ax.set_xlabel('Subclass')
            ax.set_xticks(np.arange(nsubclass))
            ax.set_xticklabels(subclasslist)
            ax.set_ylabel('Eigenvalue Ratio, $a_1/a_0$')
            ax.set_title('STAR {0}: Eigenvalue Ratios'.format(c))
            fig.savefig(outfile+'.{0}.ratio.png'.format(c))
            pylab.close(fig)
    #
    # Save output to FITS file.
    #
    if os.path.exists(outfile+'.fits'):
        os.remove(outfile+'.fits')
    hdu0 = pyfits.PrimaryHDU(fullflux)
    hdu1 = pyfits.new_table(pyfits.ColDefs([
        pyfits.Column(name='plate',format='J',array=slist['plate']),
        pyfits.Column(name='mjd',format='J',array=slist['mjd']),
        pyfits.Column(name='fiber',format='J',array=slist['fiberid']),
        pyfits.Column(name='redshift',format='D',unit='km/s',array=slist['cz'])]))
    hdulist = pyfits.HDUList([hdu0,hdu1])
    hdulist[0].header.update('OBJECT','STAR')
    hdulist[0].header.update('COEFF0',pcaflux['newloglam'][0])
    hdulist[0].header.update('COEFF1',objdloglam)
    hdulist[0].header.update('IDLUTILS','pydlutils','Version of idlutils')
    hdulist[0].header.update('SPEC2D','eigenspectra','Version of idlspec2d')
    hdulist[0].header.update('RUN2D',os.getenv('RUN2D'),'Version of 2d reduction')
    hdulist[0].header.update('RUN1D',os.getenv('RUN1D'),'Version of 1d reduction')
    for i in range(len(namearr)):
        hdulist[0].header.update("NAME%d" % i,namearr[i]+' ')
    hdulist[1].header.update('FILENAME',inputfile)
    hdulist.writeto(outfile+'.fits')
    # plot_eig(outfile+'.fits')
    return
