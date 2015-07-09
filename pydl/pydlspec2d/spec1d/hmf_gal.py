# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Create galaxy template files with HMF.
"""
from __future__ import print_function
#
def hmf_gal(**kwargs):
    """Wrapper on hmf_solve analogous to pca_gal and pca_solve.

    Parameters
    ----------
    inputfile : str, optional
        The list of spectra to use.  If not specified, $IDLSPEC2D_DIR/tempates/eigeninput_gal.dat will be used.
    wavemin : float, optional
        Minimum wavelength for the template.  If not specified 1900 Å will be used.
    wavemax : float, optional
        Maximum wavelength for the template.  If not specified 10000 Å will be used.
    K : int, optional
        Number of templates to calculate.  The default is 4.
    nonnegative : bool, optional
        If set to ``True`` use nonnegative HMF.  The default is ``False``.
    epsilon : float, optional
        Value of regularization parameter to use.  The default is 0.0, which means it is not used.
    flux : bool, optional
        If set to ``True`` make some additional QA plots of the input spectra.
        The default is ``False``.

    Returns
    -------
    None

    Notes
    -----
    Creates spEigenGal-MJD.fits and some associated QA plots.  These files
    will be created in the same directory as the inputfile (see above).

    The ``$RUN2D`` environment variable must be set.  This routine will search
    for a pickle file of the form ``eigeninput_gal_$RUN2D.dump``.  If this file
    is not found, it will be created.
    """
    import os
    import os.path
    import pickle
    import matplotlib
    matplotlib.use('Agg') # Non-interactive back-end
    import pylab
    from astropy.io import ascii, fits
    import numpy as np
    # from matplotlib.font_manager import fontManager, FontProperties
    from ...goddard.astro import get_juldate
    from ...pydlutils.image import djs_maskinterp
    from ...pydlutils.math import djs_median, find_contiguous
    from astropy.io import ascii
    from . import plot_eig, readspec, skymask, wavevector, hmf_solve
    from ..spec2d import combine1fiber
    if 'inputfile' in kwargs:
        inputfile = kwargs['inputfile']
    else:
        inputfile = os.path.join(os.getenv('IDLSPEC2D_DIR'),
            'templates','eigeninput_gal.dat')
    outdir = os.path.dirname(inputfile)
    if 'wavemin' in kwargs:
        wavemin = kwargs['wavemin']
    else:
        # Almost everything below 1900 is missing
        wavemin = 1900.0
        # wavemin = 1850.0
    if 'wavemax' in kwargs:
        wavemax = kwargs['wavemax']
    else:
        wavemax = 10000.0
    snmax = 100.0
    if 'K' in kwargs:
        K = kwargs['K']
    else:
        K = 4
    if 'nonnegative' in kwargs:
        nonnegative = kwargs['nonnegative']
    else:
        nonnegative = False
    if 'epsilon' in kwargs:
        epsilon = kwargs['epsilon']
    else:
        epsilon = 0.0
    if 'flux' in kwargs:
        plot_flux = kwargs['flux']
    else:
        plot_flux = False
    #
    # Name the output files.
    #
    jd = get_juldate()
    outfile = os.path.join(outdir,"spEigenGal-{0:d}".format(int(jd - 2400000.5)))
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
    #
    # Run combine1fiber
    #
    dump = os.path.join(outdir,'eigeninput_gal_{0}.dump'.format(os.getenv('RUN2D')))
    if os.path.exists(dump):
        print("Loading data from {0}.".format(dump))
        f = open(dump)
        foo = pickle.load(f)
        newflux = foo['newflux']
        newivar = foo['newivar']
        newloglam = foo['newloglam']
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
        nobj, npix = spplate['flux'].shape
        #
        # The redshift of each object in pixels would be logshift/objdloglam.
        #
        logshift = np.log10(1.0 + zfit)
        #
        # Determine the new wavelength mapping.
        #
        fullloglam = newloglam
        dloglam = fullloglam[1] - fullloglam[0]
        nnew = fullloglam.size
        fullflux = np.zeros((nobj,nnew),dtype='d')
        fullivar = np.zeros((nobj,nnew),dtype='d')
        #
        # Shift each spectrum to z = 0 and sample at the output wavelengths
        #
        if spplate['loglam'].ndim == 1:
            indx = spplate['loglam'] > 0
            rowloglam = loglam[indx]
        for iobj in range(nobj):
            print("OBJECT {0:5d}".format(iobj))
            if spplate['loglam'].ndim > 1:
                if spplate['loglam'].shape[0] != nobj:
                    raise ValueError('Wrong number of dimensions for loglam.')
                indx = spplate['loglam'][iobj,:] > 0
                rowloglam = spplate['loglam'][iobj,indx]
            flux1,ivar1 = combine1fiber(rowloglam-logshift[iobj],spplate['flux'][iobj,indx],
                objinvvar[iobj,indx],newloglam=fullloglam,binsz=dloglam,aesthetics='mean') # ,verbose=True)
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
        foo = dict()
        foo['newflux'] = newflux
        foo['newivar'] = newivar
        foo['newloglam'] = newloglam
        f = open(dump,'w')
        pickle.dump(foo,f)
        f.close()
    #
    # Find the columns out side of which there is no data at all
    #
    si = newflux*newivar
    zerocol = (newflux.sum(0) == 0) | (newivar.sum(0) == 0) | (si.sum(0) == 0)
    #
    # Find the largest set of contiguous pixels
    #
    goodcol = find_contiguous(~zerocol)
    newflux = newflux[:,goodcol]
    newivar = newivar[:,goodcol]
    # si = si[:,goodcol]
    newloglam = newloglam[goodcol]
    N,M = newflux.shape
    #
    # Run the HMF iteration
    #
    a,g = hmf_solve(newflux,newivar,
        K=K,nonnegative=nonnegative,epsilon=epsilon)
    #
    # Make plots
    #
    colorvec = ['k','r','g','b','m','c']
    # smallfont = FontProperties(size='xx-small');
    nspectra = N
    if plot_flux:
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
    p = ax.plot(10.0**newloglam,(newivar == 0).sum(0)/float(nspectra),'k-')
    ax.set_xlabel(r'Wavelength [$\AA$]')
    ax.set_ylabel('Fraction of spectra with missing data')
    ax.set_title('Missing Data')
    fig.savefig(outfile+'.missing.png')
    pylab.close(fig)
    # aratio10 = pcaflux['acoeff'][:,1]/pcaflux['acoeff'][:,0]
    # aratio20 = pcaflux['acoeff'][:,2]/pcaflux['acoeff'][:,0]
    # aratio30 = pcaflux['acoeff'][:,3]/pcaflux['acoeff'][:,0]
    # fig = pylab.figure(dpi=100)
    # ax = fig.add_subplot(111)
    # p = ax.plot(aratio10,aratio20,marker='None',linestyle='None')
    # for k in range(len(aratio10)):
    #     t = ax.text(aratio10[k],aratio20[k],'%04d-%04d'%(plate[k],fiber[k]),
    #         horizontalalignment='center', verticalalignment='center',
    #         color=colorvec[k%len(colorvec)],
    #         fontproperties=smallfont)
    # ax.set_xlim([aratio10.min(),aratio10.max])
    # ax.set_xlim([aratio20.min(),aratio20.max])
    # ax.set_xlabel('Eigenvalue Ratio, $a_1/a_0$')
    # ax.set_ylabel('Eigenvalue Ratio, $a_2/a_0$')
    # ax.set_title('Galaxies: Eigenvalue Ratios')
    # fig.savefig(outfile+'.a2_v_a1.png')
    # pylab.close(fig)
    # fig = pylab.figure(dpi=100)
    # ax = fig.add_subplot(111)
    # p = ax.plot(aratio20,aratio30,marker='None',linestyle='None')
    # for k in range(len(aratio10)):
    #     t = ax.text(aratio20[k],aratio30[k],'%04d-%04d'%(plate[k],fiber[k]),
    #         horizontalalignment='center', verticalalignment='center',
    #         color=colorvec[k%len(colorvec)],
    #         fontproperties=smallfont)
    # ax.set_xlim([aratio10.min(),aratio10.max])
    # ax.set_xlim([aratio20.min(),aratio20.max])
    # ax.set_xlabel('Eigenvalue Ratio, $a_2/a_0$')
    # ax.set_ylabel('Eigenvalue Ratio, $a_3/a_0$')
    # ax.set_title('Galaxies: Eigenvalue Ratios')
    # fig.savefig(outfile+'.a3_v_a2.png')
    # pylab.close(fig)
    #
    # Save output to FITS file.
    #
    if os.path.exists(outfile+'.fits'):
        os.remove(outfile+'.fits')
    hdu0 = fits.PrimaryHDU(g)
    hdu1 = fits.new_table(fits.ColDefs([
        fits.Column(name='plate',format='J',array=plate),
        fits.Column(name='mjd',format='J',array=mjd),
        fits.Column(name='fiber',format='J',array=fiber),
        fits.Column(name='redshift',format='D',array=zfit)]))
    hdulist = fits.HDUList([hdu0,hdu1])
    hdulist[0].header.update('OBJECT','GALAXY','Type of template')
    hdulist[0].header.update('COEFF0',newloglam[0],'ln(lambda) of the first spectral pixel')
    hdulist[0].header.update('COEFF1',newloglam[1]-newloglam[0],'Delta ln(lambda)')
    hdulist[0].header.update('NONNEG',nonnegative,'Was nonnegative HMF used?')
    hdulist[0].header.update('EPSILON',epsilon,'Regularization parameter used.')
    hdulist[0].header.update('IDLUTILS','pydlutils','Version of idlutils')
    hdulist[0].header.update('SPEC2D','hmf','Version of idlspec2d')
    hdulist[0].header.update('RUN2D',os.getenv('RUN2D'),'Version of 2d reduction')
    hdulist[0].header.update('RUN1D',os.getenv('RUN1D'),'Version of 1d reduction')
    # for i in range(len(pcaflux['eigenval'])):
    #     hdulist[0].header.update("EIGEN%d" % i,pcaflux['eigenval'][i])
    hdulist[1].header.update('FILENAME',inputfile,'Original input file')
    hdulist.writeto(outfile+'.fits')
    plot_eig(outfile+'.fits')
    return
#
#
#
def main():
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
    parser.add_argument('-n', '--nonnegative', action='store_true', dest='nonnegative',
        help='Use non-negative HMF method.')
    parser.add_argument('-e', '--epsilon', action='store', type=float, dest='epsilon',
        metavar='EPSILON', default=1.0, help='Set the epsilon parameter (default 1.0). Set to 0 to turn off entirely')
    parser.add_argument('-K', '--dims', action='store', type=int, dest='K',
        metavar='K', default=4, help='Set the number of functions to model (default 4).')
    #parser.add_argument('-d', '--outdir', action='store', dest='outdir',
    #    metavar='DIR', default=os.path.join(os.getenv('HOME'),'scratch'),
    #    help='Write output files to DIR.')
    parser.add_argument('-f', '--file', action='store', dest='inputfile',
        metavar='FILE', default=os.path.join(os.getenv('HOME'),'scratch','eigeninput_gal.dat'),
        help='Read input spectra and redshifts from FILE.')
    options = parser.parse_args()
    hmf_gal(inputfile=options.inputfile,dump=options.dump,K=options.K,epsilon=options.epsilon,
        nonnegative=options.nonnegative,flux=options.flux)
    return
