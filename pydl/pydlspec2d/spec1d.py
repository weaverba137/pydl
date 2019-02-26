# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""This module corresponds to the spec1d directory in idlspec2d.
"""
import glob
import os
import time
from warnings import warn
import numpy as np
from numpy.linalg import solve
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['figure.figsize'] = (16.0, 12.0)
import matplotlib.pyplot as plt
from matplotlib.font_manager import fontManager, FontProperties
from astropy import log
from astropy.io import ascii, fits
from . import Pydlspec2dException, Pydlspec2dUserWarning

#
# Used by findspec
#
findspec_cache = None


class HMF(object):
    """Class used to manage data for Heteroscedastic Matrix Factorization (HMF).

    This is a replacement for :func:`~pydl.pydlspec2d.spec1d.pca_solve`.
    It can be called with::

        hmf = HMF(spectra, invvar)
        output = hmf.solve()

    The input spectra should be pre-processed through
    :func:`~pydl.pydlspec2d.spec2d.combine1fiber`.

    Parameters
    ----------
    spectra : array-like
        The input spectral flux, assumed to have a common wavelength and
        redshift system.
    invvar : array-like
        The inverse variance of the spectral flux.
    K : :class:`int`, optional
        The number of dimensions of the factorization (default 4).
    n_iter : :class:`int`, optional
        Number of iterations.
    seed : :class:`int`, optional.
        If set, pass this value to :func:`numpy.random.seed`.
    nonnegative : :class:`bool`, optional
        Set this to ``True`` to perform nonnegative HMF.
    epsilon : :class:`float`, optional
        Regularization parameter.  Set to any non-negative float value to turn
        it on.
    verbose : :class:`bool`, optional
        If ``True``, print extra information.

    Notes
    -----
    See [1]_ and [2]_ for the original derivation of this method.

    The HMF iteration is initialized using :func:`~scipy.cluster.vq.kmeans`,
    which itself uses random numbers to initialize its state.  If you need
    to ensure reproducibility, call :func:`numpy.random.seed` before
    initializing HMF.

    The current algorithm cannot handle input data that contain *columns*
    of zeros.  Columns of this type need to be *carefully* removed from the
    input data.  This could also result in the output data having a different
    size compared to the input data.

    References
    ----------
    .. [1] `Tsalmantza, P., Decarli, R., Dotti, M., Hogg, D. W., 2011 ApJ 738, 20
        <http://adsabs.harvard.edu/abs/2011ApJ...738...20T>`_
    .. [2] `Tsalmantza, P., Hogg, D. W., 2012 ApJ 753, 122
        <http://adsabs.harvard.edu/abs/2012ApJ...753..122T>`_
    """

    def __init__(self, spectra, invvar, K=4, n_iter=None, seed=None,
                 nonnegative=False, epsilon=None, verbose=False):
        self.spectra = spectra
        self.invvar = invvar
        self.K = K
        if n_iter is None:
            if nonnegative:
                self.n_iter = 2048
            else:
                self.n_iter = 20
        else:
            self.n_iter = int(n_iter)
        self.seed = seed
        self.nonnegative = nonnegative
        self.epsilon = epsilon
        self.verbose = verbose
        self.a = None
        self.g = None
        if self.verbose:
            log.setLevel('DEBUG')
        return

    def solve(self):
        """Process the inputs.

        Returns
        -------
        :class:`dict`
            The HMF solution.
        """
        if len(self.spectra.shape) == 1:
            nobj = 1
            npix = self.spectra.shape[0]
        else:
            nobj, npix = self.spectra.shape
        log.info("Building HMF from %d object spectra.", nobj)
        fluxdict = dict()
        #
        # If there is only one object spectrum, then all we can do is return it.
        #
        if nobj == 1:
            fluxdict['flux'] = self.spectra.astype('f')
            return fluxdict
        a, g = self.iterate()
        fluxdict['acoeff'] = a
        fluxdict['flux'] = g
        return fluxdict

    def model(self):
        """Compute the model.
        """
        return np.dot(self.a, self.g)

    def resid(self):
        """Compute residuals.
        """
        return self.spectra - self.model()

    def chi(self):
        """Compute :math:`\chi`, the scaled residual.
        """
        return self.resid() * np.sqrt(self.invvar)

    def penalty(self):
        """Compute penalty for non-smoothness.
        """
        if self.epsilon is None:
            return 0.0
        return self.epsilon * np.sum(np.diff(self.g)**2)

    def badness(self):
        """Compute :math:`\chi^2`, including possible non-smoothness penalty.
        """
        return np.sum(self.chi()**2) + self.penalty()

    def normbase(self):
        """Apply standard component normalization.
        """
        return np.sqrt((self.g**2).mean(1))

    def astep(self):
        """Update for coefficients at fixed component spectra.
        """
        N, M = self.spectra.shape
        K, M = self.g.shape
        a = np.zeros((N, K), dtype=self.g.dtype)
        for i in range(N):
            Gi = np.zeros((K, K), dtype=self.g.dtype)
            for k in range(K):
                for kp in range(k, K):
                    Gi[k, kp] = np.sum(self.g[k, :] * self.g[kp, :] *
                                       self.invvar[i, :])
                    if kp > k:
                        Gi[kp, k] = Gi[k, kp]
            Fi = np.dot(self.g, self.spectra[i, :]*self.invvar[i, :])
            a[i, :] = solve(Gi, Fi)
        return a

    def gstep(self):
        """Update for component spectra at fixed coefficients.
        """
        N, M = self.spectra.shape
        N, K = self.a.shape
        g = np.zeros((K, M), dtype=self.a.dtype)
        e = np.zeros(self.g.shape, dtype=self.g.dtype)
        d = np.zeros((K, K, M), dtype=self.a.dtype)
        if self.epsilon is not None and self.epsilon > 0:
            foo = self.epsilon * np.eye(K, dtype=self.a.dtype)
            for l in range(M):
                d[:, :, l] = foo
                if l > 0 and l < M-1:
                    d[:, :, l] *= 2
            # d[:, :, 0] = foo
            # d[:, :, 1:M-1] = 2*foo
            # d[:, :, M-1] = foo
            e[:, 0] = self.epsilon*self.g[:, 1]
            e[:, 1:M-1] = self.epsilon*(self.g[:, 0:M-2] + self.g[:, 2:M])
            e[:, M-1] = self.epsilon*self.g[:, M-2]
        for j in range(M):
            Aj = np.zeros((K, K), dtype=self.a.dtype)
            for k in range(K):
                for kp in range(k, K):
                    Aj[k, kp] = np.sum(self.a[:, k] * self.a[:, kp] *
                                       self.invvar[:, j])
                    if kp > k:
                        Aj[kp, k] = Aj[k, kp]
            Aj += d[:, :, j]
            Fj = (np.dot(self.a.T, self.spectra[:, j]*self.invvar[:, j]) +
                  e[:, j])
            g[:, j] = solve(Aj, Fj)
        return g

    def astepnn(self):
        """Non-negative update for coefficients at fixed component spectra.
        """
        numerator = np.dot(self.spectra*self.invvar, self.g.T)
        denominator = np.dot(np.dot(self.a, self.g)*self.invvar, self.g.T)
        return self.a*(numerator/denominator)

    def gstepnn(self):
        """Non-negative update for component spectra at fixed coefficients.
        """
        K, M = self.g.shape
        numerator = np.dot(self.a.T, (self.spectra*self.invvar))
        if self.epsilon is not None and self.epsilon > 0:
            e = np.zeros(self.g.shape, dtype=self.g.dtype)
            e[:, 0] = self.epsilon*self.g[:, 1]
            e[:, 1:M-1] = self.epsilon*(self.g[:, 0:M-2] + self.g[:, 2:M])
            e[:, M-1] = self.epsilon*self.g[:, M-2]
            numerator += e
        denominator = np.dot(self.a.T, np.dot(self.a, self.g)*self.invvar)
        if self.epsilon is not None and self.epsilon > 0:
            d = self.epsilon*self.g.copy()
            d[:, 1:M-1] *= 2
            denominator += d
        return self.g*(numerator/denominator)

    def reorder(self):
        """Reorder and rotate basis analogous to PCA.
        """
        from numpy.linalg import eigh
        l, U = eigh(np.dot(self.a.T, self.a))
        return (np.dot(self.a, U), np.dot(U.T, self.g))

    def iterate(self):
        """Handle the HMF iteration.

        Returns
        -------
        :func:`tuple` of :class:`numpy.ndarray`
            The fitting coefficients and fitted functions, respectively.
        """
        from scipy.cluster.vq import kmeans, whiten
        from ..pydlutils.math import find_contiguous
        N, M = self.spectra.shape
        #
        # Make spectra non-negative
        #
        if self.nonnegative:
            self.spectra[self.spectra < 0] = 0
            self.invvar[self.spectra < 0] = 0
        #
        # Detect and fix very bad columns.
        #
        si = self.spectra * self.invvar
        if (self.spectra.sum(0) == 0).any():
            log.warn("Columns of zeros detected in spectra!")
        if (self.invvar.sum(0) == 0).any():
            log.warn("Columns of zeros detected in invvar!")
        if (si.sum(0) == 0).any():
            log.warn("Columns of zeros detected in spectra*invvar!")
        zerocol = ((self.spectra.sum(0) == 0) | (self.invvar.sum(0) == 0) |
                   (si.sum(0) == 0))
        n_zero = zerocol.sum()
        if n_zero > 0:
            log.warn("Found %d bad columns in input data!", n_zero)
        #
        # Find the largest set of contiguous pixels
        #
        goodcol = find_contiguous(~zerocol)
        self.spectra = self.spectra[:, goodcol]
        self.invvar = self.invvar[:, goodcol]
        # si = si[:, goodcol]
        # newloglam = fullloglam[goodcol]
        #
        # Initialize g matrix with kmeans
        #
        if self.seed is not None:
            np.random.seed(self.seed)
        whitespectra = whiten(self.spectra)
        log.debug(whitespectra[0:3, 0:3])
        self.g, foo = kmeans(whitespectra, self.K)
        self.g /= np.repeat(self.normbase(), M).reshape(self.g.shape)
        log.debug(self.g[0:3, 0:3])
        #
        # Initialize a matrix
        #
        self.a = np.outer(np.sqrt((self.spectra**2).mean(1)),
                          np.repeat(1.0/self.K, self.K))
        if self.nonnegative:
            for k in range(128):
                self.a = self.astepnn()
        #
        # Iterate!
        #
        t0 = time.time()
        for m in range(self.n_iter):
            log.info("Starting iteration #%4d.", m+1)
            if self.nonnegative:
                self.a = self.astepnn()
                self.g = self.gstepnn()
            else:
                self.a = self.astep()
                self.g = self.gstep()
                self.a, self.g = self.reorder()
            norm = self.normbase()
            self.g /= np.repeat(norm, M).reshape(self.g.shape)
            self.a = (self.a.T*np.repeat(norm, N).reshape(self.K, N)).T
            log.debug(self.a[0:3, 0:3])
            log.debug(self.g[0:3, 0:3])
            log.debug("Chi**2 after iteration #%4d = %f.", m+1, self.badness())
            log.info("The elapsed time for iteration #%4d is %6.2f s.", m+1, time.time()-t0)
        return (self.a, self.g)


def findspec(*args, **kwargs):
    """Find SDSS/BOSS spectra that match a given RA, Dec.

    Parameters
    ----------
    ra, dec : array-like, optional
        If set, the first two positional arguments will be interpreted as
        RA, Dec.
    best : :class:`bool`, optional
        If set, return only the best match for each input RA, Dec.
    infile : :class:`str`, optional
        If set, read RA, Dec data from this file.
    outfile : :class:`str`, optional
        If set, print match data to this file.
    print : :class:`bool`, optional
        If set, print the match data to the console.
    run1d : :class:`str`, optional
        Override the value of :envvar:`RUN1D`.
    run2d : :class:`str`, optional
        Override the value of :envvar:`RUN2D`.
    sdss : :class:`bool`, optional
        If set, search for SDSS-I/II spectra instead of BOSS spectra.
    searchrad : :class:`float`, optional
        Search for spectra in this radius around given RA, Dec.
        Default is 3 arcsec.
    topdir : :class:`str`, optional
        If set, override the value of :envvar:`SPECTRO_REDUX`
        or :envvar:`BOSS_SPECTRO_REDUX`.

    Returns
    -------
    :class:`dict`
        A dictionary containing plate, MJD, fiber, etc.
    """
    from .. import uniq
    from ..pydlutils.misc import struct_print
    from ..pydlutils.spheregroup import spherematch
    global findspec_cache
    #
    # Set up default values
    #
    if 'sdss' in kwargs:
        if 'topdir' in kwargs:
            topdir = kwargs['topdir']
        else:
            topdir = os.environ['SPECTRO_REDUX']
        if 'run2d' in kwargs:
            run2d = str(kwargs['run2d'])
        else:
            run2d = '26'
        run1d = ''
    else:
        if 'topdir' in kwargs:
            topdir = kwargs['topdir']
        else:
            topdir = os.environ['BOSS_SPECTRO_REDUX']
        if 'run2d' in kwargs:
            run2d = str(kwargs['run2d'])
        else:
            run2d = os.environ['RUN2D']
        if 'run1d' in kwargs:
            run1d = str(kwargs['run1d'])
        else:
            run1d = os.environ['RUN1D']
    if findspec_cache is None:
        findspec_cache = {'lasttopdir': topdir, 'plist': None}
    if (findspec_cache['plist'] is None or
        topdir != findspec_cache['lasttopdir']):
        findspec_cache['lasttopdir'] = topdir
        platelist_file = os.path.join(topdir, "platelist.fits")
        plates_files = glob.glob(os.path.join(topdir, "plates-*.fits"))
        plist = None
        if os.path.exists(platelist_file):
            platelist = fits.open(platelist_file)
            plist = platelist[1].data
            platelist.close()
        if len(plates_files) > 0:
            plates = fits.open(plates_files[0])
            plist = plates[1].data
            plates.close()
        if plist is None:
            raise Pydlspec2dException("Plate list (platelist.fits or plates-*.fits) not found in {0}.".format(topdir))
        else:
            findspec_cache['plist'] = plist
    qdone = plist.field('STATUS1D') == 'Done'
    qdone2d = plist.field('RUN2D').strip() == run2d
    if run1d == '':
        qdone1d = np.ones(plist.size, dtype='bool')
    else:
        qdone1d = plist.field('RUN1D').strip() == run1d
    qfinal = qdone & qdone2d & qdone1d
    if not qfinal.any():
        warn("No reduced plates!", Pydlspec2dUserWarning)
        return None
    idone = np.arange(plist.size)[qfinal]
    #
    # If there are positional arguments, interpret these as RA, Dec
    #
    if len(args) == 2:
        ra = args[0]
        dec = args[1]
    #
    # Read RA, Dec from infile if set
    #
    if 'infile' in kwargs:
        infile_data = ascii.read(kwargs['infile'], names=['ra', 'dec'])
        ra = infile_data["ra"].data
        dec = infile_data["dec"].data
    if 'searchrad' in kwargs:
        searchrad = float(kwargs['searchrad'])
    else:
        searchrad = 3.0/3600.0
    #
    # Create output structure
    #
    slist_type = np.dtype([('PLATE', 'i4'), ('MJD', 'i4'), ('FIBERID', 'i4'),
                          ('RA', 'f8'), ('DEC', 'f8'), ('MATCHRAD', 'f8')])
    #
    # Match all plates with objects
    #
    imatch1, itmp, dist12 = spherematch(ra, dec, plist[qfinal].field('RACEN'),
                                        plist[qfinal].field('DECCEN'),
                                        searchrad+1.55, maxmatch=0)
    if imatch1.size == 0:
        warn("No matching plates found.", Pydlspec2dUserWarning)
        return None
    imatch2 = idone[itmp]
    #
    # Read all relevant plates
    #
    try:
        n_total = plist.field('N_TOTAL')
    except KeyError:
        n_total = np.zeros(plist.size, dtype='i4') + 640
    iplate = imatch2[uniq(imatch2, imatch2.argsort())]
    i0 = 0
    plugmap = np.zeros(n_total[iplate].sum(), dtype=[('PLATE', 'i4'),
                       ('MJD', 'i4'), ('FIBERID', 'i4'),
                       ('RA', 'd'), ('DEC', 'd')])
    for i in range(iplate.size):
        spplate = readspec(plist[iplate[i]].field('PLATE'),
                           mjd=plist[iplate[i]].field('MJD'),
                           topdir=topdir,
                           run2d=run2d, run1d=run1d)
        index_to = i0 + np.arange(n_total[iplate[i]], dtype='i4')
        plugmap['PLATE'][index_to] = plist[iplate[i]].field('PLATE')
        plugmap['MJD'][index_to] = plist[iplate[i]].field('MJD')
        plugmap['FIBERID'][index_to] = spplate['plugmap']['FIBERID']
        plugmap['RA'][index_to] = spplate['plugmap']['RA']
        plugmap['DEC'][index_to] = spplate['plugmap']['DEC']
        i0 += n_total[iplate[i]]
    i1, i2, d12 = spherematch(ra, dec, plugmap['RA'], plugmap['DEC'],
                              searchrad, maxmatch=0)
    if i1.size == 0:
        warn('No matching objects found.', Pydlspec2dUserWarning)
        return None
    if 'best' in kwargs:
        #
        # Return only best match per object
        #
        slist = np.zeros(ra.size, dtype=slist_type)
        spplate = readspec(plugmap[i2]['PLATE'],
                           plugmap[i2]['FIBERID'],
                           mjd=plugmap[i2]['MJD'],
                           topdir=topdir,
                           run2d=run2d, run1d=run1d)
        sn = spplate['zans']['SN_MEDIAN']
        isort = (i1 + np.where(sn > 0, sn, 0)/(sn+1.0).max()).argsort()
        i1 = i1[isort]
        i2 = i2[isort]
        d12 = d12[isort]
        iuniq = uniq(i1)
        slist[i1[iuniq]]['PLATE'] = plugmap[i2[iuniq]]['PLATE']
        slist[i1[iuniq]]['MJD'] = plugmap[i2[iuniq]]['MJD']
        slist[i1[iuniq]]['FIBERID'] = plugmap[i2[iuniq]]['FIBERID']
        slist[i1[iuniq]]['RA'] = plugmap[i2[iuniq]]['RA']
        slist[i1[iuniq]]['DEC'] = plugmap[i2[iuniq]]['DEC']
        slist[i1[iuniq]]['MATCHRAD'] = d12[iuniq]
    else:
        #
        # Return all matches
        #
        slist = np.zeros(i1.size, dtype=slist_type)
        slist['PLATE'] = plugmap[i2]['PLATE']
        slist['MJD'] = plugmap[i2]['MJD']
        slist['FIBERID'] = plugmap[i2]['FIBERID']
        slist['RA'] = plugmap[i2]['RA']
        slist['DEC'] = plugmap[i2]['DEC']
        slist['MATCHRAD'] = d12
    #
    # Print to terminal or output file
    #
    if 'print' in kwargs:
        foo = struct_print(slist)
    if 'outfile' in kwargs:
        foo = struct_print(slist, filename=outfile)
    return slist


def latest_mjd(plate, **kwargs):
    """Find the most recent MJD associated with a plate.

    Parameters
    ----------
    plate : :class:`int` or :class:`numpy.ndarray`
        The plate(s) to examine.

    Returns
    -------
    :class:`numpy.ndarray`
        An array of MJD values for each plate.
    """
    import re
    if isinstance(plate, (int,)) or plate.shape == ():
        platevec = np.array([plate], dtype='i4')
    else:
        platevec = plate
    mjd = np.zeros(len(platevec), dtype='i4')
    mjdre = re.compile(r'spPlate-[0-9]{4}-([0-9]{5}).fits')
    unique_plates = np.unique(platevec)
    paths = spec_path(unique_plates, **kwargs)
    for p, q in zip(paths, unique_plates):
        plateglob = "{0}/spPlate-{1:04d}-*.fits".format(p, q)
        bigmjd = 0
        for f in glob.glob(plateglob):
            thismjd = int(mjdre.search(f).groups()[0])
            if thismjd > bigmjd:
                bigmjd = thismjd
        mjd[platevec == q] = bigmjd
    return mjd


def number_of_fibers(plate, **kwargs):
    """Returns the total number of fibers per plate.

    Parameters
    ----------
    plate : :class:`int` or :class:`numpy.ndarray`
        The plate(s) to examine.

    Returns
    -------
    :class:`numpy.ndarray`
        The number of fibers on each plate.
    """
    #
    # Get mjd values
    #
    if isinstance(plate, (int,)) or plate.shape == ():
        platevec = np.array([plate], dtype='i4')
    else:
        platevec = plate
    mjd = latest_mjd(plate, **kwargs)
    nfiber = np.zeros(mjd.size, dtype='i4')
    #
    # SDSS-I,II plates
    #
    nfiber[mjd < 55025] = 640
    #
    # Short circuit if we're done.
    #
    if (nfiber == 640).all():
        return nfiber
    #
    # Not all BOSS plates have 1000 fibers
    #
    if 'path' in kwargs:
        platelistpath = os.path.join(kwargs['path'], 'platelist.fits')
    else:
        platelistpath = os.path.join(os.environ['BOSS_SPECTRO_REDUX'], 'platelist.fits')
    platelist = fits.open(platelistpath)
    platentotal = platelist[1].data.field('N_TOTAL')
    plateplate = platelist[1].data.field('PLATE')
    platemjd = platelist[1].data.field('MJD')
    platerun2d = platelist[1].data.field('RUN2D')
    platerun1d = platelist[1].data.field('RUN1D')
    platelist.close()
    if 'run2d' in kwargs:
        run2d = kwargs['run2d']
    else:
        run2d = os.environ['RUN2D']
    if 'run1d' in kwargs:
        run1d = kwargs['run1d']
    else:
        run1d = os.environ['RUN1D']
    for k in range(mjd.size):
        nfiber[k] = platentotal[(plateplate == platevec[k]) &
                                (platemjd == mjd[k]) &
                                (platerun2d == run2d) &
                                (platerun1d == run1d)]
    return nfiber


def pca_solve(newflux, newivar, maxiter=0, niter=10, nkeep=3,
              nreturn=None, verbose=False):
    """Replacement for idlspec2d pca_solve.pro.

    Parameters
    ----------
    newflux : array-like
        The input spectral flux, assumed to have a common wavelength and
        redshift system.
    newivar : array-like
        The inverse variance of the spectral flux.
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
    from .. import pcomp
    from ..pydlutils.math import computechi2, djs_reject
    if verbose:
        log.setLevel('DEBUG')
    if nreturn is None:
        nreturn = nkeep
    if len(newflux.shape) == 1:
        nobj = 1
        npix = newflux.shape[0]
    else:
        nobj, npix = newflux.shape
    log.info("Building PCA from %d object spectra.", nobj)
    nzi = newivar.nonzero()
    first_nonzero = (np.arange(nobj, dtype=nzi[0].dtype),
                     np.array([nzi[1][nzi[0] == k].min() for k in range(nobj)]))
    #
    # Construct the synthetic weight vector, to be used when replacing the
    # low-S/N object pixels with the reconstructions.
    #
    synwvec = np.ones((npix,), dtype='d')
    for ipix in range(npix):
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
            flux0 = np.tile(filtflux[first_nonzero], npix).reshape(npix, nobj).transpose()
            # flux0 = np.tile(filtflux, npix).reshape(npix, nobj).transpose()
            totflux = np.absolute(filtflux - flux0).sum(1)
            goodobj = totflux > 0
            if goodobj.all():
                tmp = pcomp(filtflux.T)  # , standardize=True)
                pres = tmp.derived
                eigenval = tmp.eigenvalues
            else:
                tmp = pcomp(filtflux[goodobj, :].T)  # , standardize=True)
                pres = np.zeros((nobj, npix), dtype='d')
                pres[goodobj, :] = tmp.derived
                eigenval = np.zeros((nobj,), dtype='d')
                eigenval[goodobj] = tmp.eigenvalues
            maskivar = newivar * outmask
            sqivar = np.sqrt(maskivar)
            for iobj in range(nobj):
                out = computechi2(newflux[iobj, :], sqivar[iobj, :],
                                  pres[:, 0:nkeep])
                filtflux[iobj, :] = (maskivar[iobj, :] * newflux[iobj, :] +
                                     synwvec*out.yfit) / (maskivar[iobj, :] +
                                                          synwvec)
                acoeff[iobj, :] = out.acoeff
            log.info("The elapsed time for iteration #%2d is %6.2f s.", ipiter+1, time.time()-t0)
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


def plot_eig(filename, title='Unknown'):
    """Plot spectra from an eigenspectra/template file.

    Parameters
    ----------
    filename : :class:`str`
        Name of a FITS file containing eigenspectra/templates.
    title : :class:`str`, optional
        Title to put on the plot.

    Raises
    ------
    :exc:`ValueError`
        If an unknown template type was input in `filename`.
    """
    #
    # Set title based on filename
    #
    if title == 'Unknown':
        if filename.find('Gal') > 0:
            title = 'Galaxies: Eigenspectra'
        elif filename.find('QSO') > 0:
            title = 'QSOs: Eigenspectra'
        elif filename.find('Star') > 0:
            title = 'Stars: Eigenspectra'
        elif filename.find('CVstar') > 0:
            title = 'CV Stars: Eigenspectra'
        else:
            raise ValueError('Unknown template type!')
    base, ext = filename.split('.')
    spectrum = fits.open(filename)
    newloglam0 = spectrum[0].header['COEFF0']
    objdloglam = spectrum[0].header['COEFF1']
    spectro_data = spectrum[0].data
    spectrum.close()
    (neig, ndata) = spectro_data.shape
    newloglam = np.arange(ndata) * objdloglam + newloglam0
    lam = 10.0**newloglam
    fig = plt.figure(dpi=100)
    ax = fig.add_subplot(111)
    colorvec = ['k', 'r', 'g', 'b', 'm', 'c']
    for l in range(neig):
        p = ax.plot(lam, spectro_data[l, :],
                    colorvec[l % len(colorvec)]+'-', linewidth=1)
    ax.set_xlabel(r'Wavelength [$\AA$]')
    ax.set_ylabel('Flux [Arbitrary Units]')
    ax.set_title(title)
    # ax.set_xlim([3500.0,10000.0])
    # ax.set_ylim([-400.0,500.0])
    # fig.savefig(base+'.zoom.png')
    fig.savefig(base+'.png')
    plt.close(fig)
    return


def readspec(platein, mjd=None, fiber=None, **kwargs):
    """Read SDSS/BOSS spec2d & spec1d files.

    Parameters
    ----------
    platein : :class:`int` or :class:`numpy.ndarray`
        Plate number(s).
    mjd : :class:`int` or :class:`numpy.ndarray`, optional
        MJD numbers.  If not provided, they will be calculated by
        :func:`latest_mjd`.
    fiber : array-like, optional
        Fibers to read.  If not set, all fibers from all plates will be
        returned.
    topdir : :class:`str`, optional
        Override the value of :envvar:`BOSS_SPECTRO_REDUX`.
    run2d : :class:`str`, optional
        Override the value of :envvar:`RUN2D`.
    run1d : :class:`str`, optional
        Override the value of :envvar:`RUN1D`.
    path : :class:`str`, optional
        Override all path information with this directory name.
    align : :class:`bool`, optional
        If set, align all the spectra in wavelength.
    znum : :class:`int`, optional
        If set, return the znum-th best fit reshift fit, instead of the best.

    Returns
    -------
    :class:`dict`
        A dictionary containing the data read.
    """
    try:
        nplate = len(platein)
        plate = platein
    except TypeError:
        nplate = 1
        plate = np.array([platein], dtype='i4')
    if 'run2d' in kwargs:
        run2d = kwargs['run2d']
    else:
        run2d = os.environ['RUN2D']
    if 'run1d' in kwargs:
        run1d = kwargs['run1d']
    else:
        run1d = os.environ['RUN1D']
    if fiber is None:
        #
        # Read all fibers
        #
        nfibers = number_of_fibers(plate, **kwargs)
        total_fibers = nfibers.sum()
        platevec = np.zeros(total_fibers, dtype='i4')
        fibervec = np.zeros(total_fibers, dtype='i4')
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
            platevec = np.array(plate, dtype='i4')
        else:
            platevec = np.zeros(nfiber, dtype='i4') + plate
        if nfiber > 1:
            fibervec = np.array(fiber, dtype='i4')
        else:
            fibervec = np.zeros(nplate, dtype='i4') + fiber
    if 'mjd' is None:
        mjdvec = latest_mjd(platevec, **kwargs)
    else:
        try:
            nmjd = len(mjd)
        except TypeError:
            nmjd = 1
        if nmjd != nplate:
            raise TypeError("Plate & MJD must have the same length!")
        mjdvec = np.zeros(nplate, dtype='i4') + mjd
    #
    # Now select unique plate-mjd combinations & read them
    #
    pmjd = ((np.array(platevec, dtype='u8') << 16) +
            np.array(mjdvec, dtype='u8'))
    # log.debug(pmjd)
    upmjd = np.unique(pmjd)
    zupmjd = list(zip(upmjd >> 16, upmjd & ((1 << 16) - 1)))
    # log.debug(zupmjd)
    spplate_data = dict()
    hdunames = ('flux', 'invvar', 'andmask', 'ormask', 'disp', 'plugmap',
                'sky', 'loglam',)
    for thisplate, thismjd in zupmjd:
        # thisplate = int(p>>16)
        # thismjd = int(np.bitwise_and(p, (1<<16)-1))
        pmjdindex = ((platevec == thisplate) &
                     (mjdvec == thismjd)).nonzero()[0]
        thisfiber = fibervec[pmjdindex]
        # log.debug(type(thisplate), type(thismjd))
        # log.debug(repr(thisfiber))
        # log.debug(type(thisfiber))
        pmjdstr = "{0:04d}-{1:05d}".format(int(thisplate), int(thismjd))
        if 'path' in kwargs:
            sppath = [kwargs['path']]
        else:
            sppath = spec_path(thisplate, run2d=run2d)
        spfile = os.path.join(sppath[0], "spPlate-{0}.fits".format(pmjdstr))
        log.info(spfile)
        spplate = fits.open(spfile)
        #
        # Get wavelength coefficients from primary header
        #
        npix = spplate[0].header['NAXIS1']
        c0 = spplate[0].header['COEFF0']
        c1 = spplate[0].header['COEFF1']
        coeff0 = np.zeros(thisfiber.size, dtype='d') + c0
        coeff1 = np.zeros(thisfiber.size, dtype='d') + c1
        loglam0 = c0 + c1*np.arange(npix, dtype='d')
        loglam = np.resize(loglam0, (thisfiber.size, npix))
        #
        # Read the data images
        #
        for k in range(len(hdunames)):
            if hdunames[k] == 'loglam':
                tmp = loglam
            else:
                try:
                    tmp = spplate[k].data[thisfiber-1, :]
                except IndexError:
                    tmp = spplate[k].data[thisfiber-1]
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
                    for c in spplate[k].columns.names:
                        spplate_data['plugmap'][c] = tmp[c]
                else:
                    spplate_data[hdunames[k]] = tmp
            else:
                #
                # Append data
                #
                if k == 0:
                    allpmjdindex = np.concatenate((allpmjdindex, pmjdindex))
                    if 'align' in kwargs:
                        mincoeff0 = min(allcoeff0)
                        if mincoeff0 == 0 and coeff0[0] > 0:
                            allcoeff0 = coeff0[0]
                            allcoeff1 = coeff1[1]
                        if mincoeff0 > 0 and coeff0[0] == 0:
                            coeff0 = mincoeff0
                            coeff1 = allcoeff1[0]
                        ps = np.floor((coeff0[0] - mincoeff0)/coeff1[0] + 0.5)
                        if ps > 0:
                            coeff0 = coeff0 - ps*coeff1
                        else:
                            allcoeff0 = allcoeff0 + ps*allcoeff1
                    else:
                        ps = 0
                    allcoeff0 = np.concatenate((allcoeff0, coeff0))
                    allcoeff1 = np.concatenate((allcoeff1, coeff1))
                if hdunames[k] == 'plugmap':
                    for c in spplate[5].columns.names:
                        spplate_data['plugmap'][c] = np.concatenate(
                            (spplate_data['plugmap'][c], tmp[c]))
                else:
                    spplate_data[hdunames[k]] = spec_append(spplate_data[hdunames[k]], tmp, pixshift=ps)
        spplate.close()
        #
        # Read photoPlate information, if available
        #
        photofile = os.path.join(sppath[0],
                                 "photoPlate-{0}.fits".format(pmjdstr))
        if not os.path.exists(photofile):
            #
            # Hmm, maybe this is an SDSS-I,II plate
            #
            photofile = os.path.join(os.environ['SPECTRO_MATCH'], run2d,
                                     os.path.basename(os.environ['PHOTO_RESOLVE']),
                                     "{0:04d}".format(int(thisplate)),
                                     "photoPlate-{0}.fits".format(pmjdstr))
        if os.path.exists(photofile):
            photop = fits.open(photofile)
            tmp = photop[1].data[thisfiber-1]
            if 'tsobj' not in spplate_data:
                spplate_data['tsobj'] = dict()
                for c in photop[1].columns.names:
                    spplate_data['tsobj'][c] = tmp[c]
            else:
                for c in photop[1].columns.names:
                    spplate_data['tsobj'][c] = np.concatenate(
                        (spplate_data['tsobj'][c], tmp[c]))
            photop.close()

        #
        # Read redshift information, if available.
        #
        if 'znum' in kwargs:
            zfile = os.path.join(sppath[0], run1d,
                                 "spZall-{0}.fits".format(pmjdstr))
        else:
            zfile = os.path.join(sppath[0], run1d,
                                 "spZbest-{0}.fits".format(pmjdstr))
        if os.path.exists(zfile):
            spz = fits.open(zfile)
            if 'znum' in kwargs:
                nper = spz[0].header['DIMS0']
                zfiber = (thisfiber-1)*nper + kwargs['znum'] - 1
            else:
                zfiber = thisfiber
            tmp = spz[1].data[zfiber-1]
            if 'zans' not in spplate_data:
                spplate_data['zans'] = dict()
                for c in spz[1].columns.names:
                    spplate_data['zans'][c] = tmp[c]
            else:
                for c in spz[1].columns.names:
                    spplate_data['zans'][c] = np.concatenate(
                        (spplate_data['zans'][c], tmp[c]))
            spz.close()
    #
    # Reorder the data.  At this point allpmjdindex is an index for which
    # fiber[allpmjdindex] == spplate['plugmap']['FIBERID'], so we have to
    # reverse this mapping.
    #
    j = allpmjdindex.argsort()
    for k in spplate_data:
        if isinstance(spplate_data[k], dict):
            for c in spplate_data[k]:
                if spplate_data[k][c].ndim == 2:
                    spplate_data[k][c] = spplate_data[k][c][j, :]
                else:
                    spplate_data[k][c] = spplate_data[k][c][j]
        else:
            spplate_data[k] = spplate_data[k][j, :]
    allcoeff0 = allcoeff0[j]
    allcoeff1 = allcoeff1[j]
    #
    # If necessary, recompute the wavelengths
    #
    nfibers, npixmax = spplate_data['flux'].shape
    if 'align' in kwargs:
        loglam0 = allcoeff0[0] + allcoeff1[1]*np.arange(npixmax, dtype='d')
        spplate_data['loglam'] = np.resize(loglam0, (nfibers, npixmax))
    return spplate_data


def skymask(invvar, andmask, ormask=None, ngrow=2):
    """Mask regions where sky-subtraction errors are expected to dominate.

    Parameters
    ----------
    invvar : :class:`numpy.ndarray`
        Inverse variance.
    andmask : :class:`numpy.ndarray`
        An "and" mask.  For historical reasons, this input is ignored.
    ormask : :class:`numpy.ndarray`, optional
        An "or" mask.  Although technically this is optional, if it is
        not supplied, this function will have no effect.
    ngrow : :class:`int`, optional
        Expand bad areas by this number of pixels.

    Returns
    -------
    :class:`numpy.ndarray`
        The `invvar` multiplied by the bad areas.
    """
    from ..pydlutils.sdss import sdss_flagval
    from .. import smooth
    nrows, npix = invvar.shape
    badmask = np.zeros(invvar.shape, dtype='i4')
    badskychi = sdss_flagval('SPPIXMASK', 'BADSKYCHI')
    redmonster = sdss_flagval('SPPIXMASK', 'REDMONSTER')
    # brightsky = sdss_flagval('SPPIXMASK', 'BRIGHTSKY')
    if ormask is not None:
        badmask = badmask | ((ormask & badskychi) != 0)
        badmask = badmask | ((ormask & redmonster) != 0)
        # badmask = badmask | ((andmask & brightsky) != 0)
    if ngrow > 0:
        width = 2*ngrow + 1
        for k in range(nrows):
            badmask[k, :] = smooth(badmask[k, :]*width, width, True) > 0
    return invvar * (1 - badmask)


def spec_append(spec1, spec2, pixshift=0):
    """Append the array spec2 to the array spec1 & return a new array.

    If the dimension of these arrays is the same, then append as [spec1,spec2].
    If not, increase the size of the smaller array & fill with zeros.

    Parameters
    ----------
    spec1, spec2 : :class:`numpy.ndarray`
        Append `spec2` to `spec1`.
    pixshift : :class:`int`, optional
        If `pixshift` is set to a positive integer, `spec2` will be padded with
        `pixshift` zeros on the left side.  If `pixshift` is set to a
        negative integer, `spec1` will be padded with ``abs(pixshift)`` zeros
        on the left side.  If not set, all zeros will be padded on the right
        side.

    Returns
    -------
    :class:`numpy.ndarray`
        A new array containing both `spec1` and `spec2`.
    """
    nrows1, npix1 = spec1.shape
    nrows2, npix2 = spec2.shape
    nrows = nrows1+nrows2
    nadd1 = 0
    nadd2 = 0
    if pixshift != 0:
        if pixshift < 0:
            nadd1 = -pixshift
        else:
            nadd2 = pixshift
    maxpix = max(npix1 + nadd1, npix2 + nadd2)
    spec3 = np.zeros((nrows, maxpix), dtype=spec1.dtype)
    spec3[0:nrows1, nadd1:nadd1+npix1] = spec1
    spec3[nrows1:nrows, nadd2:nadd2+npix2] = spec2
    return spec3


def spec_path(plate, path=None, topdir=None, run2d=None):
    """Return the directory containing spPlate files.

    Parameters
    ----------
    plate : :class:`int` or :class:`numpy.ndarray`
        The plate(s) to examine.
    path : :class:`str`, optional
        If set, `path` becomes the full path for every plate. In other words,
        it completely short-circuits this function.
    topdir : :class:`str`, optional
        Used to override the value of :envvar:`BOSS_SPECTRO_REDUX`.
    run2d : :class:`str`, optional
        Used to override the value of :envvar:`RUN2D`.

    Returns
    -------
    :class:`list`
        A list of directories, one for each plate.

    Raises
    ------
    :exc:`KeyError`
        If environment variables are not supplied.
    """
    if isinstance(plate, (int,)) or plate.shape == ():
        platevec = np.array([plate], dtype='i4')
    else:
        platevec = plate
    if path is None:
        if run2d is None:
            run2d = os.environ['RUN2D']
        if topdir is None:
            env = "SPECTRO_REDUX"
            try:
                ir = int(run2d)
            except ValueError:
                env = 'BOSS_SPECTRO_REDUX'
            topdir = os.environ[env]
    paths = list()
    for p in platevec:
        if path is not None:
            paths.append(path)
        else:
            paths.append(os.path.join(topdir, run2d, '{0:04d}'.format(p)))
    return paths


def preprocess_spectra(flux, ivar, loglam=None, zfit=None, aesthetics='mean',
                       newloglam=None, wavemin=None, wavemax=None,
                       verbose=False):
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
    verbose : :class:`bool`, optional
        If ``True``, print extra information.

    Returns
    -------
    :func:`tuple` of :class:`numpy.ndarray`
        The resampled flux, inverse variance and wavelength solution,
        respectively.
    """
    from .spec2d import combine1fiber
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
            log.info("OBJECT %5d", iobj)
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
        return (fullflux, fullivar, fullloglam)


def template_metadata(inputfile, verbose=False):
    """Read template metadata from file.

    Parameters
    ----------
    inputfile : :class:`str`
        Name of a Parameter file containing the input data and metadata.
    verbose : :class:`bool`, optional
        If ``True``, print lots of extra information.

    Returns
    -------
    :func:`tuple`
        A tuple containing the list of input spectra and a dictionary
        containing other metadata.
    """
    from ..pydlutils.yanny import yanny
    if verbose:
        log.setLevel('DEBUG')
    if not os.path.exists(inputfile):
        raise Pydlspec2dException("Could not find {0}!".format(inputfile))
    log.debug("Reading input data from %s.", inputfile)
    par = yanny(inputfile)
    required_metadata = {'object': str, 'method': str, 'aesthetics': str,
                         'run2d': str, 'run1d': str,
                         'wavemin': float, 'wavemax': float, 'snmax': float,
                         'niter': int, 'nkeep': int, 'minuse': int}
    metadata = dict()
    for key in required_metadata:
        try:
            metadata[key] = required_metadata[key](par[key])
            log.debug('%s = %s', key, par[key])
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
    if metadata['method'].lower() == 'hmf':
        required_hmf_metadata = {'nonnegative': lambda x: bool(int(x)),
                                 'epsilon': float}
        for key in required_hmf_metadata:
            try:
                metadata[key] = required_hmf_metadata[key](par[key])
            except KeyError:
                raise KeyError('The {0} keyword was not found in {1}!'.format(key, inputfile))
            except ValueError:
                raise ValueError('The {0} keyword has invalid value, {0}!'.format(key, par[key]))
    return (slist, metadata)


def template_input(inputfile, dumpfile, flux=False, verbose=False):
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
    """
    import pickle
    from astropy.constants import c as cspeed
    from .. import uniq
    from .. import __version__ as pydl_version
    from ..goddard.astro import get_juldate
    from ..pydlutils.image import djs_maskinterp
    from ..pydlutils.math import djs_median
    #
    # Logging
    #
    if verbose:
        log.setLevel('DEBUG')
    #
    # Read metadata.
    #
    slist, metadata = template_metadata(inputfile)
    #
    # Name the output files.
    #
    jd = get_juldate()
    outfile = "spEigen{0}-{1:d}".format(metadata['object'].title(), int(jd - 2400000.5))
    #
    # Read the input spectra
    #
    if os.path.exists(dumpfile):
        log.info("Loading data from %s.", dumpfile)
        with open(dumpfile) as f:
            inputflux = pickle.load(f)
        newflux = inputflux['newflux']
        newivar = inputflux['newivar']
        newloglam = inputflux['newloglam']
    else:
        if metadata['object'].lower() == 'star':
            spplate = readspec(slist.plate, mjd=slist.mjd, fiber=slist.fiberid,
                               align=True)
        else:
            spplate = readspec(slist.plate, mjd=slist.mjd, fiber=slist.fiberid)
        #
        # Insist that all of the requested spectra exist.
        #
        missing = spplate['plugmap']['FIBERID'] == 0
        if missing.any():
            imissing = missing.nonzero()[0]
            for k in imissing:
                log.error("Missing plate=%d mjd=%d fiberid=%d",
                          slist.plate[k], slist.mjd[k], slist.fiberid[k])
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
        if metadata['object'].lower() == 'star':
            newloglam = spplate['loglam'][0, :]
        else:
            newloglam = wavevector(np.log10(metadata['wavemin']),
                                   np.log10(metadata['wavemax']), binsz=objdloglam)
        try:
            zfit = slist.zfit
        except AttributeError:
            zfit = slist.cz/cspeed.to('km / s').value
        #
        # Shift to common wavelength grid.
        #
        newflux, newivar, newloglam = preprocess_spectra(spplate['flux'],
                                                         objinvvar,
                                                         loglam=spplate['loglam'],
                                                         zfit=zfit,
                                                         newloglam=newloglam,
                                                         aesthetics=metadata['aesthetics'],
                                                         verbose=verbose)
        #
        # Dump input fluxes to a file for debugging purposes.
        #
        if not os.path.exists(dumpfile):
            with open(dumpfile, 'w') as f:
                inputflux = {'newflux': newflux, 'newivar': newivar,
                             'newloglam': newloglam}
                pickle.dump(inputflux, f)
    #
    # Solve.
    #
    if metadata['object'].lower() == 'qso':
        pcaflux = template_qso(metadata, newflux, newivar, verbose)
    elif metadata['object'].lower() == 'star':
        pcaflux = template_star(metadata, newloglam, newflux, newivar,
                                slist, outfile, verbose)
    else:
        if metadata['method'].lower() == 'pca':
            pcaflux = pca_solve(newflux, newivar,
                                niter=metadata['niter'],
                                nkeep=metadata['nkeep'],
                                verbose=verbose)
        elif metadata['method'].lower() == 'hmf':
            hmf = HMF(newflux, newivar,
                      K=metadata['nkeep'],
                      n_iter=metadata['niter'],
                      nonnegative=metadata['nonnegative'],
                      epsilon=metadata['epsilon'],
                      verbose=verbose)
            pcaflux = hmf.solve()
        else:
            raise ValueError("Unknown method: {0}!".format(metadata['method']))
    pcaflux['newflux'] = newflux
    pcaflux['newivar'] = newivar
    pcaflux['newloglam'] = newloglam
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
    # Make plots
    #
    colorvec = ['k', 'r', 'g', 'b', 'm', 'c']
    smallfont = FontProperties(size='xx-small')
    nspectra = pcaflux['newflux'].shape[0]
    #
    # Plot input spectra
    #
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
    #
    # Missing data diagnostic.
    #
    fig = plt.figure(dpi=100)
    ax = fig.add_subplot(111)
    p = ax.plot(10.0**pcaflux['newloglam'], (pcaflux['newivar'] == 0).sum(0)/float(nspectra), 'k-')
    ax.set_xlabel(r'Wavelength [$\AA$]')
    ax.set_ylabel('Fraction of spectra with missing data')
    ax.set_title('Missing Data')
    ax.grid(True)
    fig.savefig(outfile+'.missing.png')
    plt.close(fig)
    #
    # usemask diagnostic
    #
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
    #
    # This type of figure isn't really meaningful for stars.
    #
    if metadata['object'].lower() != 'star':
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
    if not pydl_version:
        pydl_version = 'git'
    hdu0.header['OBJECT'] = (objtypes[metadata['object']], 'Type of template')
    hdu0.header['COEFF0'] = (pcaflux['newloglam'][0], 'Wavelength zeropoint')
    hdu0.header['COEFF1'] = (pcaflux['newloglam'][1]-pcaflux['newloglam'][0], 'Delta wavelength')
    hdu0.header['IDLUTILS'] = ('pydl-{0}'.format(pydl_version), 'Version of idlutils')
    hdu0.header['SPEC2D'] = ('pydl-{0}'.format(pydl_version), 'Version of idlspec2d')
    hdu0.header['RUN2D'] = (os.environ['RUN2D'], 'Version of 2d reduction')
    hdu0.header['RUN1D'] = (os.environ['RUN1D'], 'Version of 1d reduction')
    hdu0.header['FILENAME'] = (inputfile, 'Input file')
    hdu0.header['METHOD'] = (metadata['method'].upper(), 'Method used')
    if metadata['method'].lower() == 'hmf':
        hdu0.header['NONNEG'] = (metadata['nonnegative'], 'Was nonnegative HMF used?')
        hdu0.header['EPSILON'] = (metadata['epsilon'], 'Regularization parameter used.')
    # for i in range(len(namearr)):
    #     hdu0.header["NAME{0:d}".format(i)] = namearr[i]+' '
    c = [fits.Column(name='plate', format='J', array=slist.plate),
         fits.Column(name='mjd', format='J', array=slist.mjd),
         fits.Column(name='fiberid', format='J', array=slist.fiberid)]
    if metadata['object'].lower() == 'star':
        c.append(fits.Column(name='cz', format='D', unit='km/s',
                 array=slist.cz))
        for i, name in enumerate(pcaflux['namearr']):
            hdu0.header['NAME{0:d}'.format(i)] = (name, 'Name of class {0:d}.'.format(i))
    else:
        c.append(fits.Column(name='zfit', format='D', array=slist.zfit))
    hdu1 = fits.BinTableHDU.from_columns(fits.ColDefs(c))
    hdulist = fits.HDUList([hdu0, hdu1])
    hdulist.writeto(outfile+'.fits')
    if metadata['object'].lower() != 'star':
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


def template_qso(metadata, newflux, newivar, verbose=False):
    """Run PCA or HMF on QSO spectra.

    Historically, QSO templates were comptuted one at a time instead of
    all at once.

    Parameters
    ----------
    metadata : :class:`dict`
        Dictionary containing metadata about the spectra.
    newflux : :class:`~numpy.ndarray`
        Flux shifted onto common wavelength.
    newivar : :class:`~numpy.ndarray`
        Inverse variances of the fluxes.
    verbose : :class:`bool`, optional
        If ``True``, print lots of extra information.

    Returns
    -------
    :class:`dict`
        A dictonary containing flux, eigenvalues, etc.
    """
    from ..pydlutils.math import computechi2
    if metadata['object'].lower() != 'qso':
        raise Pydlspec2dException("You appear to be passing the wrong kind of object to template_qso()!")
    if len(newflux.shape) == 1:
        nobj = 1
        npix = newflux.shape[0]
    else:
        nobj, npix = newflux.shape
    objflux = newflux.copy()
    for ikeep in range(metadata['nkeep']):
        log.info("Solving for eigencomponent #%d of %d", ikeep+1, nkeep)
        if metadata['method'].lower() == 'pca':
            pcaflux1 = pca_solve(objflux, newivar,
                                 niter=metadata['niter'], nkeep=1,
                                 verbose=verbose)
        elif metadata['method'].lower() == 'hmf':
            hmf = HMF(objflux, newivar,
                      K=metadata['nkeep'],
                      n_iter=metadata['niter'],
                      nonnegative=metadata['nonnegative'],
                      epsilon=metadata['epsilon'],
                      verbose=verbose)
            pcaflux1 = hmf.solve()
        else:
            raise ValueError("Unknown method: {0}!".format(metadata['method']))
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
            # objloglam = None
            if ikeep == 0:
                objflux = newflux - np.outer(acoeff, pcaflux['flux'])
            else:
                objflux = newflux - np.dot(acoeff, pcaflux['flux'])
            # objflux = newflux - np.outer(acoeff,pcaflux['flux'])
            # objinvvar = pcaflux1['newivar']
            pcaflux['acoeff'] = acoeff
    return pcaflux


def template_star(metadata, newloglam, newflux, newivar, slist, outfile,
                  verbose=False):
    """Run PCA or HMF on stellar spectra of various classes.

    Parameters
    ----------
    metadata : :class:`dict`
        Dictionary containing metadata about the spectra.
    newloglam : :class:`~numpy.ndarray`
        The wavelength array, used only for plots.
    newflux : :class:`~numpy.ndarray`
        Flux shifted onto common wavelength.
    newivar : :class:`~numpy.ndarray`
        Inverse variances of the fluxes.
    slist : :class:`~numpy.recarray`
        The list of objects, containing stellar class information.
    outfile : :class:`str`
        The base name of output file, used for plots.
    verbose : :class:`bool`, optional
        If ``True``, print lots of extra information.

    Returns
    -------
    :class:`dict`
        A dictonary containing flux, eigenvalues, etc.
    """
    from .. import uniq
    from ..pydlutils.image import djs_maskinterp
    if metadata['object'].lower() != 'star':
        raise Pydlspec2dException("You appear to be passing the wrong kind of object to template_star()!")
    #
    # Find the list of unique star types
    #
    isort = np.argsort(slist['class'])
    classlist = slist['class'][isort[uniq(slist['class'][isort])]]
    #
    # Loop over each star type
    #
    npix, nstars = newflux.shape
    pcaflux = dict()
    pcaflux['namearr'] = list()
    for c in classlist:
        #
        # Find the subclasses for this stellar type
        #
        log.info("Finding eigenspectra for Stellar class %s.", c)
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
        if metadata['method'].lower() == 'pca':
            pcaflux1 = pca_solve(newflux[indx, :], newivar[indx, :],
                                niter=metadata['niter'], nkeep=nkeep,
                                verbose=verbose)
        elif metadata['method'].lower() == 'hmf':
            hmf = HMF(newflux[indx, :], newivar[indx, :],
                      K=metadata['nkeep'],
                      n_iter=metadata['niter'],
                      nonnegative=metadata['nonnegative'],
                      epsilon=metadata['epsilon'],
                      verbose=verbose)
            pcaflux1 = hmf.solve()
        else:
            raise ValueError("Unknown method: {0}!".format(metadata['method']))
        #
        # Some star templates are generated from only one spectrum,
        # and these will not have a usemask set.
        #
        if 'usemask' not in pcaflux1:
            pcaflux1['usemask'] = np.zeros((npix,), dtype='i4') + nindx
        #
        # Interpolate over bad flux values in the middle of a spectrum,
        # and set fluxes to zero at the blue+red ends of the spectrum
        #
        # minuse = 1 # ?
        minuse = np.floor((nindx+1) / 3.0)
        qbad = pcaflux1['usemask'] < minuse
        #
        # Interpolate over all bad pixels
        #
        for j in range(nkeep):
            pcaflux1['flux'][j, :] = djs_maskinterp(pcaflux1['flux'][j, :],
                                                    qbad, const=True)
        #
        # Set bad pixels at the very start or end of the spectrum to zero
        # instead.
        #
        npix = qbad.size
        igood = (~qbad).nonzero()[0]
        if qbad[0]:
            pcaflux1['flux'][:, 0:igood[0]-1] = 0
        if qbad[npix-1]:
            pcaflux1['flux'][:, igood[::-1][0]+1:npix] = 0
        #
        # Re-normalize the first eigenspectrum to a mean of 1
        #
        norm = pcaflux1['flux'][0, :].mean()
        pcaflux1['flux'] /= norm
        if 'acoeff' in pcaflux1:
            pcaflux1['acoeff'] *= norm
        #
        # Now loop through each stellar subclass and reconstruct
        # an eigenspectrum for that subclass
        #
        thesesubclassnum = np.zeros(thesesubclass.size, dtype='i4')
        colorvec = ['k', 'r', 'g', 'b', 'm', 'c']
        smallfont = FontProperties(size='xx-small')
        fig = plt.figure(dpi=100)
        ax = fig.add_subplot(111)
        for isub in range(nsubclass):
            ii = (thesesubclass == subclasslist[isub]).nonzero()[0]
            thesesubclassnum[ii] = isub
            if nkeep == 1:
                thisflux = pcaflux1['flux'][0, :]
            else:
                aratio = pcaflux1['acoeff'][ii, 1]/pcaflux1['acoeff'][ii, 0]
                #
                # np.median(foo) is equivalent to MEDIAN(foo,/EVEN)
                #
                thisratio = np.median(aratio)
                thisflux = (pcaflux1['flux'][0, :] +
                            thisratio.astype('f') * pcaflux1['flux'][1, :])
            #
            # Plot spectra
            #
            plotflux = thisflux/thisflux.max()
            ax.plot(10.0**newloglam, plotflux,
                    "{0}-".format(colorvec[isub % len(colorvec)]),
                    linewidth=1)
            if isub == 0:
                ax.set_xlabel(r'Wavelength [$\AA$]')
                ax.set_ylabel('Flux [arbitrary units]')
                ax.set_title('STAR {0}: Eigenspectra Reconstructions'.format(c))
            t = ax.text(10.0**newloglam[-1], plotflux[-1],
                        subclasslist[isub],
                        horizontalalignment='right', verticalalignment='center',
                        color=colorvec[isub % len(colorvec)], fontproperties=smallfont)
            fig.savefig(outfile+'.{0}.png'.format(c))
            plt.close(fig)
            if 'flux' in pcaflux:
                pcaflux['flux'] = np.vstack((pcaflux['flux'], thisflux))
                # pcaflux['acoeff'] = np.vstack((pcaflux['acoeff'], pcaflux1['acoeff']))
                # pcaflux['usemask'] = np.vstack((pcaflux['usemask'], pcaflux1['usemask']))
            else:
                pcaflux['flux'] = thisflux
                # pcaflux['acoeff'] = pcaflux1['acoeff']
                # pcaflux['usemask'] = pcaflux1['usemask']
            pcaflux['namearr'].append(subclasslist[isub])
    return pcaflux


def template_input_main():  # pragma: no cover
    """Entry point for the compute_templates script.

    Returns
    -------
    :class:`int`
        An integer suitable for passing to :func:`sys.exit`.
    """
    #
    # Imports for main()
    #
    import sys
    from argparse import ArgumentParser

    # Get home directory in platform-independent way
    home_dir = os.path.expanduser('~')
    #
    # Get Options
    #
    parser = ArgumentParser(description="Compute spectral templates.",
                            prog=os.path.basename(sys.argv[0]))
    parser.add_argument('-d', '--dump', action='store', dest='dump',
                        metavar='FILE',
                        default=os.path.join(home_dir, 'scratch', 'templates', 'compute_templates.dump'),
                        help='Dump data to a pickle file (default: %(default)s).')
    parser.add_argument('-F', '--flux', action='store_true', dest='flux',
                        help='Plot input spectra.')
    parser.add_argument('-f', '--file', action='store', dest='inputfile',
                        metavar='FILE',
                        default=os.path.join(home_dir, 'scratch', 'templates', 'compute_templates.par'),
                        help='Read input spectra and redshifts from FILE (default: %(default)s).')
    parser.add_argument('-v', '--verbose', action='store_true', dest='verbose',
                        help='Print lots of extra information.')
    options = parser.parse_args()
    template_input(options.inputfile, options.dump, options.flux, options.verbose)
    return 0


def wavevector(minfullwave, maxfullwave, zeropoint=3.5, binsz=1.0e-4,
               wavemin=None):
    """Return an array of wavelengths.

    Parameters
    ----------
    minfullwave : :class:`float`
        Minimum wavelength.
    maxfullwave : :class:`float`
        Maximum wavelength.
    zeropoint : :class:`float`, optional
        Offset of the input wavelength values.
    binsz : :class:`float`, optional
        Separation between wavelength values.
    wavemin : :class:`float`, optional
        If this is set the values of `minfullwave` and `zeropoint` are ignored.

    Returns
    -------
    :class:`numpy.ndarray`
        Depending on the values of `minfullwave`, `binsz`, etc., the resulting
        array could be interpreted as an array of wavelengths or an array of
        log(wavelength).
    """
    if wavemin is not None:
        spotmin = 0
        spotmax = int((maxfullwave - wavemin)/binsz)
        wavemax = spotmax * binsz + wavemin
    else:
        spotmin = int((minfullwave - zeropoint)/binsz) + 1
        spotmax = int((maxfullwave - zeropoint)/binsz)
        wavemin = spotmin * binsz + zeropoint
        wavemax = spotmax * binsz + zeropoint
    nfinalpix = spotmax - spotmin + 1
    finalwave = np.arange(nfinalpix, dtype='d') * binsz + wavemin
    return finalwave
