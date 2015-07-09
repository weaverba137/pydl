# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Create galaxy template files.
"""
from __future__ import print_function
#
def pca_gal(**kwargs):
    """Wrapper on pca_solve to handle galaxy eigenspectra.

    Parameters
    ----------
    inputfile : str, optional
        The list of spectra to use.  If not specified, $IDLSPEC2D_DIR/tempates/eigeninput_gal.dat will be used.
    wavemin : float, optional
        Minimum wavelength for the template.  If not specified 1850 Å will be used.
    wavemax : float, optional
        Maximum wavelength for the template.  If not specified 10000 Å will be used.
    niter : int, optional
        Number of iterations.  The default is 10.
    dump : str, optional
        If set, save input data in a Python pickle file.
    flux : bool, optional
        If set to ``True`` make some additional QA plots of the input spectra.
        The default is ``False``.

    Returns
    -------
    None

    Notes
    -----
    Creates spEigenGal-MJD.fits and some associated QA plots.
    """
    import os
    import os.path
    import pickle
    import matplotlib
    matplotlib.use('Agg') # Non-interactive back-end
    import pylab
    from astropy.io import ascii, fits
    import numpy as np
    from matplotlib.font_manager import fontManager, FontProperties
    from ...goddard.astro import get_juldate
    from ...pydlutils.image import djs_maskinterp
    from ...pydlutils.math import djs_median
    from . import pca_solve, plot_eig, readspec, skymask, wavevector
    if 'inputfile' in kwargs:
        inputfile = kwargs['inputfile']
    else:
        inputfile = os.path.join(os.getenv('IDLSPEC2D_DIR'),
            'templates','eigeninput_gal.dat')
    if 'wavemin' in kwargs:
        wavemin = kwargs['wavemin']
    else:
        wavemin = 1850.0
    if 'wavemax' in kwargs:
        wavemax = kwargs['wavemax']
    else:
        wavemax = 10000.0
    snmax = 100.0
    if 'niter' in kwargs:
        niter = kwargs['niter']
    else:
        niter = 10
    nkeep = 4
    minuse = 10
    #
    # Name the output files.
    #
    jd = get_juldate()
    outfile = "spEigenGal-%d" % int(jd - 2400000.5)
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
    if 'dump' in kwargs:
        dumpfile = kwargs['dump']
    else:
        dumpfile = 'this-file-does-not-exist'
    if os.path.exists(dumpfile):
        print("Loading data from {0}.".format(dumpfile))
        f = open(dumpfile)
        pcaflux = pickle.load(f)
        f.close()
    else:
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
            ax.set_title('Galaxies: Input Spectra %4d-%4d' % (istart+1,iend+1))
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
    ax.set_title('Galaxies: Eigenvalue Ratios')
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
    ax.set_title('Galaxies: Eigenvalue Ratios')
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
    hdulist[0].header.update('OBJECT','GALAXY')
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
#
#
#
def pca_gal_main(): # pragma: no cover
    """Gets called when the script is executed.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Side-Effects
    ------------
    Creates a FITS file and some PNG plots.
    """
    #
    # Imports for main()
    #
    import os
    import os.path
    import sys
    from astropy.utils.compat import argparse
    #
    # Get Options
    #
    parser = argparse.ArgumentParser(description=__doc__,prog=os.path.basename(sys.argv[0]))
    parser.add_argument('-F', '--flux', action='store_true', dest='flux',
        help='Plot input spectra.')
    parser.add_argument('-d', '--dump', action='store', dest='dump',
        metavar='FILE', default=os.path.join(os.getenv('HOME'),'scratch','eigeninput_gal.dump'),
        help='Dump data to a pickle file.')
    #parser.add_argument('-n', '--nonnegative', action='store_true', dest='nonnegative',
    #    help='Use non-negative HMF method.')
    #parser.add_argument('-e', '--epsilon', action='store', type=float, dest='epsilon',
    #    metavar='EPSILON', default=1.0, help='Set the epsilon parameter (default 1.0). Set to 0 to turn off entirely')
    #parser.add_argument('-K', '--dims', action='store', type=int, dest='K',
    #    metavar='K', default=4, help='Set the number of functions to model (default 4).')
    #parser.add_argument('-d', '--outdir', action='store', dest='outdir',
    #    metavar='DIR', default=os.path.join(os.getenv('HOME'),'scratch'),
    #    help='Write output files to DIR.')
    parser.add_argument('-f', '--file', action='store', dest='inputfile',
        metavar='FILE', default=os.path.join(os.getenv('HOME'),'scratch','eigeninput_gal.dat'),
        help='Read input spectra and redshifts from FILE.')
    options = parser.parse_args()
    pca_gal(inputfile=options.inputfile,dump=options.dump,
        flux=options.flux)
    return 0
