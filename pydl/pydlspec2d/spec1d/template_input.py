# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-


def template_input(options):
    """Collect spectra and pass them to PCA or HMF solvers to compute
    spectral templates.

    Parameters
    ----------
    options : parsed options
        A parsed set of command-line options.

    Returns
    -------
    None
    """
    import os
    import pickle
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams['figure.figsize'] = (16.0, 12.0)
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import fontManager, FontProperties
    from astropy.io import fits
    from astropy import log
    from warnings import warn
    from . import pca_solve, plot_eig, readspec, skymask, wavevector
    from .. import Pydlspec2dException, Pydlspec2dUserWarning
    from ... import uniq
    from ... import __version__ as pydl_version
    from ...goddard.astro import get_juldate
    from ...pydlutils.image import djs_maskinterp
    from ...pydlutils.math import djs_median
    from ...pydlutils.yanny import yanny
    #
    # Read input data
    #
    par = yanny(options.inputfile)
    try:
        objtype = par['object']
    except KeyError:
        raise KeyError('The object keyword was not found in {0}!'.format(options.inputfile))
    try:
        method = par['method']
    except KeyError:
        raise KeyError('The method keyword was not found in {0}!'.format(options.inputfile))
    try:
        wavemin = float(par['wavemin'])
    except KeyError:
        raise KeyError('The wavemin keyword was not found in {0}!'.format(options.inputfile))
    except ValueError:
        raise ValueError('The wavemin keyword has invalid value, {0}!'.format(par['wavemin']))
    try:
        wavemax = float(par['wavemax'])
    except KeyError:
        raise KeyError('The wavemax keyword was not found in {0}!'.format(options.inputfile))
    except ValueError:
        raise ValueError('The wavemax keyword has invalid value, {0}!'.format(par['wavemax']))
    try:
        snmax = float(par['snmax'])
    except KeyError:
        raise KeyError('The snmax keyword was not found in {0}!'.format(options.inputfile))
    except ValueError:
        raise ValueError('The snmax keyword has invalid value, {0}!'.format(par['snmax']))
    try:
        niter = int(par['niter'])
    except KeyError:
        raise KeyError('The niter keyword was not found in {0}!'.format(options.inputfile))
    except ValueError:
        raise ValueError('The niter keyword has invalid value, {0}!'.format(par['niter']))
    try:
        nkeep = int(par['nkeep'])
    except KeyError:
        raise KeyError('The nkeep keyword was not found in {0}!'.format(options.inputfile))
    except ValueError:
        raise ValueError('The nkeep keyword has invalid value, {0}!'.format(par['nkeep']))
    try:
        minuse = int(par['minuse'])
    except KeyError:
        raise KeyError('The minuse keyword was not found in {0}!'.format(options.inputfile))
    except ValueError:
        raise ValueError('The minuse keyword has invalid value, {0}!'.format(par['minuse']))
    cspeed = 2.99792458e5
    #
    # Name the output files.
    #
    jd = get_juldate()
    outfile = "spEigen{0}-{1:d}".format(objtype.title(), int(jd - 2400000.5))
    #
    # Read the input spectra
    #
    if os.path.exists(options.dump):
        log.info("Loading data from {0}.".format(options.dump))
        with open(options.dump) as f:
            pcaflux = pickle.load(f)
    else:
        slist = par['EIGENOBJ']
        spplate = readspec(slist.plate, mjd=slist.mjd, fiber=slist.fiberid)
        #
        # Insist that all of the requested spectra exist.
        #
        missing = spplate['plugmap']['FIBERID'] == 0
        if missing.any():
            imissing = missing.nonzero()[0]
            for k in imissing:
                log.error("Missing plate={0:d} mjd={1:d} fiberid={2:d}".format(
                        slist.plate[k], slist.mjd[k], slist.fiberid[k]))
            raise ValueError("{0:d} missing object(s).".format(missing.sum()))
        #
        # Do not fit where the spectrum may be dominated by sky-sub residuals.
        #
        objinvvar = skymask(spplate['invvar'], spplate['andmask'],
                            spplate['ormask'])
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
        try:
            objdloglam = float(par['binsz'])
        except:
            objdloglam = spplate['loglam'][0, 1] - spplate['loglam'][0, 0]
        newloglam = wavevector(np.log10(wavemin), np.log10(wavemax),
            binsz=objdloglam)
        #
        # Do PCA solution.
        #
        pcaflux = pca_solve(spplate['flux'], objinvvar, spplate['loglam'],
            slist.zfit, niter=niter, nkeep=nkeep, newloglam=newloglam,
            aesthetics='mean')
        #
        # Fill in bad data with a running median of the good data.
        #
        qgood = pcaflux['usemask'] >= minuse
        medflux = np.zeros(pcaflux['flux'].shape, dtype=pcaflux['flux'].dtype)
        if not qgood.all():
            for i in range(nkeep):
                medflux[i, qgood] = djs_median(pcaflux['flux'][i, qgood],
                    width=51, boundary='nearest')
                medflux[i, :] = djs_maskinterp(medflux[i, :], ~qgood, const=True)
            pcaflux['flux'][:, ~qgood] = medflux[:, ~qgood]
        #
        # Dump input fluxes to a file for debugging purposes.
        #
        if not os.path.exists(options.dump):
            with open(kwargs['dump'], 'w') as f:
                pickle.dump(pcaflux, f)
    #
    # Make plots
    #
    colorvec = ['k', 'r', 'g', 'b', 'm', 'c']
    smallfont = FontProperties(size='xx-small')
    nspectra = pcaflux['newflux'].shape[0]
    if options.flux:
        nfluxes = 30
        separation = 5.0
        nplots = nspectra/nfluxes
        if nspectra % nfluxes > 0:
            nplots += 1
        for k in range(nplots):
            istart = k*nfluxes
            iend = min(istart+nfluxes, nspectra) - 1
            fig = plt.figure(dpi=100)
            ax = fig.add_subplot(111)
            for l in range(istart, iend+1):
                p = ax.plot(10.0**pcaflux['newloglam'],
                    pcaflux['newflux'][l, :] + separation*(l % nfluxes),
                    colorvec[l % len(colorvec)]+'-',
                    linewidth=1)
            ax.set_xlabel(r'Wavelength [$\AA$]')
            ax.set_ylabel(r'Flux [$\mathsf{10^{-17} erg\, cm^{-2} s^{-1} \AA^{-1}}$] + Constant')
            ax.set_title('Input Spectra {0:04d}-{1:04d}'.format(istart+1, iend+1))
            ax.set_ylim(pcaflux['newflux'][istart, :].min(), pcaflux['newflux'][iend-1, :].max()+separation*(nfluxes-1))
            fig.savefig('{0}.flux.{1:04d}-{2:04d}.png'.format(outfile, istart+1, iend+1))
            plt.close(fig)
    fig = plt.figure(dpi=100)
    ax = fig.add_subplot(111)
    p = ax.plot(10.0**pcaflux['newloglam'], (pcaflux['newivar'] == 0).sum(0)/float(nspectra), 'k-')
    ax.set_xlabel(r'Wavelength [$\AA$]')
    ax.set_ylabel('Fraction of spectra with missing data')
    ax.set_title('Missing Data')
    fig.savefig(outfile+'.missing.png')
    plt.close(fig)
    aratio10 = pcaflux['acoeff'][:, 1]/pcaflux['acoeff'][:, 0]
    aratio20 = pcaflux['acoeff'][:, 2]/pcaflux['acoeff'][:, 0]
    aratio30 = pcaflux['acoeff'][:, 3]/pcaflux['acoeff'][:, 0]
    fig = plt.figure(dpi=100)
    ax = fig.add_subplot(111)
    p = ax.plot(aratio10, aratio20, marker='None', linestyle='None')
    for k in range(len(aratio10)):
        t = ax.text(aratio10[k], aratio20[k],
            '{0:04d}-{1:04d}'.format(slist.plate[k], slist.fiberid[k]),
            horizontalalignment='center', verticalalignment='center',
            color=colorvec[k % len(colorvec)],
            fontproperties=smallfont)
    # ax.set_xlim([aratio10.min(), aratio10.max])
    # ax.set_xlim([aratio20.min(), aratio20.max])
    ax.set_xlabel('Eigenvalue Ratio, $a_1/a_0$')
    ax.set_ylabel('Eigenvalue Ratio, $a_2/a_0$')
    ax.set_title('Eigenvalue Ratios')
    fig.savefig(outfile+'.a2_v_a1.png')
    plt.close(fig)
    fig = plt.figure(dpi=100)
    ax = fig.add_subplot(111)
    p = ax.plot(aratio20, aratio30, marker='None', linestyle='None')
    for k in range(len(aratio10)):
        t = ax.text(aratio20[k], aratio30[k],
            '{0:04d}-{1:04d}'.format(slist.plate[k], slist.fiberid[k]),
            horizontalalignment='center', verticalalignment='center',
            color=colorvec[k % len(colorvec)],
            fontproperties=smallfont)
    # ax.set_xlim([aratio10.min(), aratio10.max])
    # ax.set_xlim([aratio20.min(), aratio20.max])
    ax.set_xlabel('Eigenvalue Ratio, $a_2/a_0$')
    ax.set_ylabel('Eigenvalue Ratio, $a_3/a_0$')
    ax.set_title('Eigenvalue Ratios')
    fig.savefig(outfile+'.a3_v_a2.png')
    plt.close(fig)
    #
    # Save output to FITS file.
    #
    if os.path.exists(outfile+'.fits'):
        os.remove(outfile+'.fits')
    hdu0 = fits.PrimaryHDU(fullflux)
    objtypes = {'gal': 'GALAXY', 'qso': 'QSO', 'star': 'STAR'}
    hdu0.header.update('OBJECT', objtypes['objtype'], 'Type of template')
    hdu0.header.update('COEFF0', pcaflux['newloglam'][0])
    hdu0.header.update('COEFF1', objdloglam)
    hdu0.header.update('IDLUTILS', 'pydl-{0}'.format(pydl_version), 'Version of idlutils')
    hdu0.header.update('SPEC2D', 'pydl-{0}'.format(pydl_version), 'Version of idlspec2d')
    hdu0.header.update('RUN2D', os.getenv('RUN2D'),
                            'Version of 2d reduction')
    hdu0.header.update('RUN1D', os.getenv('RUN1D'),
                            'Version of 1d reduction')
    hdu0.header.update('FILENAME', inputfile)
    hdu0.header.update('METHOD', method.upper(), 'Method used')
    if method == 'hmf':
        hdu0.header.update('NONNEG', nonnegative, 'Was nonnegative HMF used?')
        hdu0.header.update('EPSILON', epsilon, 'Regularization parameter used.')
    # for i in range(len(namearr)):
    #     hdu0.header.update("NAME{0:d}".format(i), namearr[i]+' ')
    c = [fits.Column(name='plate', format='J', array=slist.plate),
         fits.Column(name='mjd', format='J', array=slist.mjd),
         fits.Column(name='fiberid', format='J', array=slist.fiberid)]
    if objtype == 'star':
        c.append(fits.Column(name='cz', format='D', unit='km/s',
            array=slist.cz))
    else:
        c.append(fits.Column(name='zfit', format='D', array=slist.zfit))
    hdu1 = fits.new_table(fits.ColDefs(c))
    hdulist = fits.HDUList([hdu0, hdu1])
    hdulist.writeto(outfile+'.fits')
    # plot_eig(outfile+'.fits')
    return


def template_input_main():  # pragma: no cover
    """Entry point for the compute_templates script.

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
    import sys
    from astropy.utils.compat import argparse
    #
    # Get Options
    #
    parser = argparse.ArgumentParser(description="Compute spectral templates.",
        prog=os.path.basename(sys.argv[0]))
    # parser.add_argument('-D', '--dims', action='store', type=int, dest='K',
    #     metavar='K', default=4, help='Set the number of functions to model (default 4).')
    parser.add_argument('-d', '--dump', action='store', dest='dump',
        metavar='FILE', default=os.path.join(os.getenv('HOME'), 'scratch', 'compute_templates.dump'),
        help='Dump data to a pickle file.')
    # parser.add_argument('-e', '--epsilon', action='store', type=float, dest='epsilon',
    #     metavar='EPSILON', default=1.0, help='Set the epsilon parameter (default 1.0). Set to 0 to turn off entirely')
    parser.add_argument('-F', '--flux', action='store_true', dest='flux',
        help='Plot input spectra.')
    parser.add_argument('-f', '--file', action='store', dest='inputfile',
        metavar='FILE', default=os.path.join(os.getenv('HOME'), 'scratch', 'compute_templates.par'),
        help='Read input spectra and redshifts from FILE.')
    # parser.add_argument('-n', '--nonnegative', action='store_true', dest='nonnegative',
    #     help='Use non-negative HMF method.')
    options = parser.parse_args()
    template_input(options)
    return 0


if __name__ == "__main__":
    from sys import exit
    exit(template_input_main())
