# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-


def preprocess_spectra(flux, ivar, loglam=None, zfit=None, aesthetics='mean',
                       newloglam=None, wavemin=None, wavemax=None,
                       good_columns=False, verbose=False):
    """Handle the processing of input spectra through the
    :func:`~pydl.pydlspec2d.spec2d.combine1fiber` stage.

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
    good_columns : :class:`bool`, optional
        If ``True``, ensure that the data contain no columns that are
        all zero.
    verbose : :class:`bool`, optional
        If ``True``, print extra information.

    Returns
    -------
    :func:`tuple` of :class:`numpy.ndarray`
        The resampled flux, inverse variance and wavelength solution,
        respectively.
    """
    import numpy as np
    from astropy import log
    from . import wavevector
    from ..spec2d import combine1fiber
    from ...pydlutils.math import find_contiguous
    if verbose:
        log.setLevel('DEBUG')
    if len(flux.shape) == 1:
        nobj = 1
        npix = flux.shape[0]
    else:
        nobj, npix = flux.shape
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
        if newloglam is None:
            raise ValueError("newloglam must be set if loglam is not!")
        return (flux, ivar, newloglam)
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
                                         flux[iobj, indx], fullloglam,
                                         objivar=ivar[iobj, indx],
                                         binsz=dloglam, aesthetics=aesthetics,
                                         verbose=verbose)
            fullflux[iobj, :] = flux1
            fullivar[iobj, :] = ivar1
        #
        # Find the columns out side of which there is no data at all.
        #
        if good_columns:
            si = fullflux*fullivar
            zerocol = ((newflux.sum(0) == 0) | (newivar.sum(0) == 0) |
                       (si.sum(0) == 0))
            #
            # Find the largest set of contiguous pixels
            #
            goodcol = find_contiguous(~zerocol)
            newflux = fullflux[:, goodcol]
            newivar = fullivar[:, goodcol]
            # si = si[:, goodcol]
            newloglam = fullloglam[goodcol]
            return (newflux, newivar, newloglam)
        else:
            return (fullflux, fullivar, fullloglam)


def template_input(inputfile, dumpfile, flux, verbose):
    """Collect spectra and pass them to PCA or HMF solvers to compute
    spectral templates.

    This function replaces the various ``PCA_GAL()``, ``PCA_STAR()``, etc.,
    functions from idlspec2d.

    Parameters
    ----------
    inputfile : :class:`str`
        Name of a Parameter file containing the input data and metadata.
    dumpfile : :class:`str`
        Name of a Pickle file used to store intermediate data.
    flux : :class:`bool`, optional
        If ``True``, plot the individual input spectra.
    verbose : :class:`bool`, optional
        If ``True``, print lots of extra information.

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
    from astropy.constants import c as cspeed
    from astropy import log
    from warnings import warn
    from . import hmf_solve, pca_solve, plot_eig, readspec, skymask, wavevector
    from .. import Pydlspec2dException, Pydlspec2dUserWarning
    from ... import uniq
    from ... import __version__ as pydl_version
    from ...goddard.astro import get_juldate
    from ...pydlutils.image import djs_maskinterp
    from ...pydlutils.math import djs_median
    from ...pydlutils.yanny import yanny
    #
    # Logging
    #
    if verbose:
        log.setLevel('DEBUG')
    #
    # Read input data
    #
    if not os.path.exists(inputfile):
        raise Pydlspec2dException("Could not find {0}!".format(inputfile))
    log.debug("Reading input data from {0}.".format(inputfile))
    par = yanny(inputfile)
    required_metadata = {'object': str, 'method': str, 'aesthetics': str,
                         'run2d': str, 'run1d': str,
                         'wavemin': float, 'wavemax': float, 'snmax': float,
                         'niter': int, 'nkeep': int, 'minuse': int}
    metadata = dict()
    for key in required_metadata:
        try:
            metadata[key] = required_metadata[key](par[key])
            log.debug('{0} = {1}'.format(key, par[key]))
        except KeyError:
            raise KeyError('The {0} keyword was not found in {1}!'.format(key, inputfile))
        except ValueError:
            raise ValueError('The {0} keyword has invalid value, {0}!'.format(key, par[key]))
    slist = par['EIGENOBJ']
    for r in ('run2d', 'run1d'):
        try:
            metadata['orig_'+r] = os.environ[r.upper()]
        except KeyError:
            metadata['orig_'+r] = None
        os.environ[r.upper()] = metadata[r]
    good_columns = False
    if metadata['method'].lower() == 'hmf':
        good_columns = True
        nonnegative = bool(metadata['nonnegative'])
        epsilon = float(metadata['epsilon'])
    #
    # Name the output files.
    #
    jd = get_juldate()
    outfile = "spEigen{0}-{1:d}".format(metadata['object'].title(), int(jd - 2400000.5))
    #
    # Read the input spectra
    #
    if os.path.exists(dumpfile):
        log.info("Loading data from {0}.".format(dumpfile))
        with open(dumpfile) as f:
            pcaflux = pickle.load(f)
    else:
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
        ifix = spplate['flux']**2 * objinvvar > metadata['snmax']**2
        if ifix.any():
            objinvvar[ifix.nonzero()] = (metadata['snmax']/spplate['flux'][ifix.nonzero()])**2
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
        newloglam = wavevector(np.log10(metadata['wavemin']),
                               np.log10(metadata['wavemax']), binsz=objdloglam)
        #
        # Do PCA solution.
        #
        newflux, newivar, newloglam = preprocess_spectra(spplate['flux'],
                                                         objinvvar,
                                                         loglam=spplate['loglam'],
                                                         zfit=slist.zfit,
                                                         newloglam=newloglam,
                                                         aesthetics=metadata['aesthetics'],
                                                         good_columns=good_columns,
                                                         verbose=verbose)
        if metadata['method'].lower() == 'pca':
            if metadata['object'].lower() == 'qso':
                #
                # Solve for one component at a time, like the old pca_qso.pro
                #
                objflux = newflux.copy()
                for ikeep in range(nkeep):
                    log.info("Solving for eigencomponent #{0:d} of {1:d}".format(ikeep+1, nkeep))
                    pcaflux1 = pca_solve(objflux, newivar,
                                         niter=metadata['niter'], nkeep=1,
                                         verbose=verbose)
                    if ikeep == 0:
                        #
                        # Create new pcaflux dict
                        #
                        pcaflux = dict()
                        for k in pcaflux1:
                            pcaflux[k] = pcaflux1[k].copy()
                    else:
                        #
                        # Add to existing dict
                        #
                        # for k in pcaflux1:
                        #     pcaflux[k] = np.vstack((pcaflux[k],pcaflux1[k]))
                        pcaflux['flux'] = np.vstack((pcaflux['flux'], pcaflux1['flux']))
                        pcaflux['eigenval'] = np.concatenate((pcaflux['eigenval'], pcaflux1['eigenval']))
                    #
                    # Re-solve for the coefficients using all PCA components so far
                    #
                    acoeff = np.zeros((nobj, ikeep+1), dtype=pcaflux1['acoeff'].dtype)
                    for iobj in range(nobj):
                        out = computechi2(newflux[iobj, :],
                                          np.sqrt(pcaflux1['newivar'][iobj, :]),
                                          pcaflux['flux'].T)
                        acoeff[iobj, :] = out['acoeff']
                    #
                    # Prevent re-binning of spectra on subsequent calls to pca_solve()
                    #
                    objloglam = None
                    if ikeep == 0:
                        objflux = newflux - np.outer(acoeff, pcaflux['flux'])
                    else:
                        objflux = newflux - np.dot(acoeff, pcaflux['flux'])
                    # objflux = newflux - np.outer(acoeff,pcaflux['flux'])
                    # objinvvar = pcaflux1['newivar']
                    pcaflux['acoeff'] = acoeff
            else:
                #
                # Do a normal simultaneous PCA solution
                #
                pcaflux = pca_solve(newflux, newivar,
                                    niter=metadata['niter'],
                                    nkeep=metadata['nkeep'],
                                    verbose=verbose)
        elif metadata['method'].lower() == 'hmf':
            pcaflux = hmf_solve(newflux, newivar,
                                K=metadata['nkeep'],
                                nonnegative=nonnegative, epsilon=epsilon,
                                verbose=verbose)
        else:
            raise ValueError("Unknown method: {0}!".format(metadata['method']))
        #
        # Fill in bad data with a running median of the good data.
        # The presence of boundary='nearest' means that this code snippet
        # was never meant to be called!  In other words it should always
        # be the case that qgood.all() is True.
        #
        if 'usemask' in pcaflux:
            qgood = pcaflux['usemask'] >= metadata['minuse']
            if not qgood.all():
                warn("Would have triggered djs_median replacement!", Pydlspec2dUserWarning)
            if False:
                medflux = np.zeros(pcaflux['flux'].shape, dtype=pcaflux['flux'].dtype)
                for i in range(metadata['nkeep']):
                    medflux[i, qgood] = djs_median(pcaflux['flux'][i, qgood],
                                                   width=51, boundary='nearest')
                    medflux[i, :] = djs_maskinterp(medflux[i, :], ~qgood, const=True)
                pcaflux['flux'][:, ~qgood] = medflux[:, ~qgood]
        #
        # Dump input fluxes to a file for debugging purposes.
        #
        pcaflux['newflux'] = newflux
        pcaflux['newivar'] = newivar
        pcaflux['newloglam'] = newloglam
        if not os.path.exists(dumpfile):
            with open(dumpfile, 'w') as f:
                pickle.dump(pcaflux, f)
    #
    # Make plots
    #
    colorvec = ['k', 'r', 'g', 'b', 'm', 'c']
    smallfont = FontProperties(size='xx-small')
    nspectra = pcaflux['newflux'].shape[0]
    if flux:
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
    ax.grid(True)
    fig.savefig(outfile+'.missing.png')
    plt.close(fig)
    if 'usemask' in pcaflux:
        fig = plt.figure(dpi=100)
        ax = fig.add_subplot(111)
        p = ax.semilogy(10.0**pcaflux['newloglam'][pcaflux['usemask'] > 0],
                        pcaflux['usemask'][pcaflux['usemask'] > 0], 'k-',
                        10.0**pcaflux['newloglam'],
                        np.zeros(pcaflux['newloglam'].shape,
                        dtype=pcaflux['newloglam'].dtype) + metadata['minuse'],
                        'k--')
        ax.set_xlabel(r'Wavelength [$\AA$]')
        ax.set_ylabel('Usemask')
        ax.set_title('UseMask')
        ax.grid(True)
        fig.savefig(outfile+'.usemask.png')
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
    hdu0 = fits.PrimaryHDU(pcaflux['flux'])
    objtypes = {'gal': 'GALAXY', 'qso': 'QSO', 'star': 'STAR'}
    hdu0.header['OBJECT'] = (objtypes[metadata['object']], 'Type of template')
    hdu0.header['COEFF0'] = (pcaflux['newloglam'][0], 'Wavelength zeropoint')
    hdu0.header['COEFF1'] = (pcaflux['newloglam'][1]-pcaflux['newloglam'][0], 'Delta wavelength')
    hdu0.header['IDLUTILS'] = ('pydl-{0}'.format(pydl_version), 'Version of idlutils')
    hdu0.header['SPEC2D'] = ('pydl-{0}'.format(pydl_version), 'Version of idlspec2d')
    hdu0.header['RUN2D'] = (os.getenv('RUN2D'), 'Version of 2d reduction')
    hdu0.header['RUN1D'] = (os.getenv('RUN1D'), 'Version of 1d reduction')
    hdu0.header['FILENAME'] = (inputfile, 'Input file')
    hdu0.header['METHOD'] = (metadata['method'].upper(), 'Method used')
    if metadata['method'] == 'hmf':
        hdu0.header['NONNEG'] = (nonnegative, 'Was nonnegative HMF used?')
        hdu0.header['EPSILON'] = (epsilon, 'Regularization parameter used.')
    # for i in range(len(namearr)):
    #     hdu0.header["NAME{0:d}".format(i)] = namearr[i]+' '
    c = [fits.Column(name='plate', format='J', array=slist.plate),
         fits.Column(name='mjd', format='J', array=slist.mjd),
         fits.Column(name='fiberid', format='J', array=slist.fiberid)]
    if metadata['object'] == 'star':
        c.append(fits.Column(name='cz', format='D', unit='km/s',
                 array=slist.cz))
    else:
        c.append(fits.Column(name='zfit', format='D', array=slist.zfit))
    hdu1 = fits.BinTableHDU.from_columns(fits.ColDefs(c))
    hdulist = fits.HDUList([hdu0, hdu1])
    hdulist.writeto(outfile+'.fits')
    plot_eig(outfile+'.fits')
    #
    # Clean up
    #
    for r in ('run2d', 'run1d'):
        if metadata['orig_'+r] is None:
            del os.environ[r.upper()]
        else:
            os.environ[r.upper()] = metadata['orig_'+r]
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

    # Get home directory in platform-independent way
    home_dir = os.path.expanduser('~')
    #
    # Get Options
    #
    parser = argparse.ArgumentParser(description="Compute spectral templates.",
                                     prog=os.path.basename(sys.argv[0]))
    parser.add_argument('-d', '--dump', action='store', dest='dump',
                        metavar='FILE',
                        default=os.path.join(home_dir, 'scratch', 'templates', 'compute_templates.dump'),
                        help='Dump data to a pickle file.')
    parser.add_argument('-F', '--flux', action='store_true', dest='flux',
                        help='Plot input spectra.')
    parser.add_argument('-f', '--file', action='store', dest='inputfile',
                        metavar='FILE',
                        default=os.path.join(home_dir, 'scratch', 'templates', 'compute_templates.par'),
                        help='Read input spectra and redshifts from FILE.')
    parser.add_argument('-v', '--verbose', action='store_true', dest='verbose',
                        help='Print lots of extra information.')
    options = parser.parse_args()
    template_input(options.inputfile, options.dump, options.flux, options.verbose)
    return 0


if __name__ == "__main__":
    from sys import exit
    exit(template_input_main())
