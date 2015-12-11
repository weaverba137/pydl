# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""This module corresponds to the bspline directory in idlutils.
"""
import numpy as np


class bspline(object):
    """Bspline class.

    Functions in the bspline library are implemented as methods on this
    class.

    Parameters
    ----------
    x : :class:`numpy.ndarray`
        The data.
    nord : :class:`int`, optional
        To be documented.
    npoly : :class:`int`, optional
        To be documented.
    bkpt : :class:`numpy.ndarray`, optional
        To be documented.
    bkspread : :class:`float`, optional
        To be documented.
    verbose : :class:`bool`, optional.
        If ``True`` print extra information.

    Attributes
    ----------
    breakpoints
        To be documented.
    nord
        To be documented.
    npoly
        To be documented.
    mask
        To be documented.
    coeff
        To be documented.
    icoeff
        To be documented.
    xmin
        To be documented.
    xmax
        To be documented.
    funcname
        To be documented.
    """

    def __init__(self, x, nord=4, npoly=1, bkpt=None, bkspread=1.0,
                 verbose=False, **kwargs):
        """Init creates an object whose attributes are similar to the
        structure returned by the create_bspline function.
        """
        #
        # Set the breakpoints.
        #
        if bkpt is None:
            startx = x.min()
            rangex = x.max() - startx
            if 'placed' in kwargs:
                w = ((kwargs['placed'] >= startx) &
                     (kwargs['placed'] <= startx+rangex))
                if w.sum() < 2:
                    bkpt = np.arange(2, dtype='f') * rangex + startx
                else:
                    bkpt = kwargs['placed'][w]
            elif 'bkspace' in kwargs:
                nbkpts = int(rangex/kwargs['bkspace']) + 1
                if nbkpts < 2:
                    nbkpts = 2
                tempbkspace = rangex/float(nbkpts-1)
                bkpt = np.arange(nbkpts, dtype='f')*tempbkspace + startx
            elif 'nbkpts' in kwargs:
                nbkpts = kwargs['nbkpts']
                if nbkpts < 2:
                    nbkpts = 2
                tempbkspace = rangex/float(nbkpts-1)
                bkpt = np.arange(nbkpts, dtype='f') * tempbkspace + startx
            elif 'everyn' in kwargs:
                npkpts = max(nx/kwargs['everyn'], 1)
                if nbkpts == 1:
                    xspot = [0]
                else:
                    xspot = int(nx/(nbkpts-1)) * np.arange(nbkpts, dtype='i4')
                bkpt = x[xspot].astype('f')
            else:
                raise ValueError('No information for bkpts.')
        imin = bkpt.argmin()
        imax = bkpt.argmax()
        if x.min() < bkpt[imin]:
            if verbose:
                print('Lowest breakpoint does not cover lowest x value: changing.')
            bkpt[imin] = x.min()
        if x.max() > bkpt[imax]:
            if verbose:
                print('Highest breakpoint does not cover highest x value: changing.')
            bkpt[imax] = x.max()
        nshortbkpt = bkpt.size
        fullbkpt = bkpt.copy()
        if nshortbkpt == 1:
            bkspace = np.float32(bkspread)
        else:
            bkspace = (bkpt[1] - bkpt[0]) * np.float32(bkspread)
        for i in np.arange(1, nord, dtype=np.float32):
            fullbkpt = np.insert(fullbkpt, 0, bkpt[0]-bkspace*i)
            fullbkpt = np.insert(fullbkpt, fullbkpt.shape[0],
                                 bkpt[nshortbkpt-1] + bkspace*i)
        #
        # Set the attributes
        #
        nc = fullbkpt.size - nord
        self.breakpoints = fullbkpt
        self.nord = nord
        self.npoly = npoly
        self.mask = np.ones((fullbkpt.size,), dtype='bool')
        if npoly > 1:
            self.coeff = np.zeros((npoly, nc), dtype='d')
            self.icoeff = np.zeros((npoly, nc), dtype='d')
        else:
            self.coeff = np.zeros((nc,), dtype='d')
            self.icoeff = np.zeros((nc,), dtype='d')
        self.xmin = 0.0
        self.xmax = 1.0
        self.funcname = 'legendre'
        return

    def fit(self, xdata, ydata, invvar, x2=None):
        """Calculate a B-spline in the least-squares sense.

        Fit is based on two variables: x which is sorted and spans a large range
        where bkpts are required y which can be described with a low order
        polynomial.

        Parameters
        ----------
        xdata : :class:`numpy.ndarray`
            Independent variable.
        ydata : :class:`numpy.ndarray`
            Dependent variable.
        invvar : :class:`numpy.ndarray`
            Inverse variance of `ydata`.
        x2 : :class:`numpy.ndarray`, optional
            Orthogonal dependent variable for 2d fits.

        Returns
        -------
        :func:`tuple`
            A tuple containing an integer error code, and the evaluation of the
            b-spline at the input values.  An error code of -2 is a failure,
            -1 indicates dropped breakpoints, 0 is success, and positive
            integers indicate ill-conditioned breakpoints.
        """
        goodbk = self.mask[self.nord:]
        nn = goodbk.sum()
        if nn < self.nord:
            yfit = np.zeros(ydata.shape, dtype='f')
            return (-2, yfit)
        nfull = nn * self.npoly
        bw = self.npoly * self.nord
        a1, lower, upper = self.action(xdata, x2=x2)
        foo = np.tile(invvar, bw).reshape(bw, invvar.size).transpose()
        a2 = a1 * foo
        alpha = np.zeros((bw, nfull+bw), dtype='d')
        beta = np.zeros((nfull+bw,), dtype='d')
        bi = np.arange(bw, dtype='i4')
        bo = np.arange(bw, dtype='i4')
        for k in range(1, bw):
            bi = np.append(bi, np.arange(bw-k, dtype='i4')+(bw+1)*k)
            bo = np.append(bo, np.arange(bw-k, dtype='i4')+bw*k)
        for k in range(nn-self.nord+1):
            itop = k*self.npoly
            ibottom = min(itop, nfull) + bw - 1
            ict = upper[k] - lower[k] + 1
            if ict > 0:
                work = np.dot(a1[lower[k]:upper[k]+1, :].T, a2[lower[k]:upper[k]+1, :])
                wb = np.dot(ydata[lower[k]:upper[k]+1], a2[lower[k]:upper[k]+1, :])
                alpha.T.flat[bo+itop*bw] += work.flat[bi]
                beta[itop:ibottom+1] += wb
        min_influence = 1.0e-10 * invvar.sum() / nfull
        errb = cholesky_band(alpha, mininf=min_influence)  # ,verbose=True)
        if isinstance(errb[0], int) and errb[0] == -1:
            a = errb[1]
        else:
            yfit, foo = self.value(xdata, x2=x2, action=a1, upper=upper, lower=lower)
            return (self.maskpoints(errb[0]), yfit)
        errs = cholesky_solve(a, beta)
        if isinstance(errs[0], int) and errs[0] == -1:
            sol = errs[1]
        else:
            #
            # It is not possible for this to get called, because cholesky_solve
            # has only one return statement, & that statement guarantees that
            # errs[0] == -1
            #
            yfit, foo = self.value(xdata, x2=x2, action=a1, upper=upper, lower=lower)
            return (self.maskpoints(errs[0]), yfit)
        if self.npoly > 1:
            self.icoeff[:, goodbk] = np.array(a[0, 0:nfull].reshape(self.npoly, nn), dtype=a.dtype)
            self.coeff[:, goodbk] = np.array(sol[0:nfull].reshape(self.npoly, nn), dtype=sol.dtype)
        else:
            self.icoeff[goodbk] = np.array(a[0, 0:nfull], dtype=a.dtype)
            self.coeff[goodbk] = np.array(sol[0:nfull], dtype=sol.dtype)
        yfit, foo = self.value(xdata, x2=x2, action=a1, upper=upper, lower=lower)
        return (0, yfit)

    def action(self, x, x2=None):
        """Construct banded bspline matrix, with dimensions [ndata, bandwidth].

        Parameters
        ----------
        x : :class:`numpy.ndarray`
            Independent variable.
        x2 : :class:`numpy.ndarray`, optional
            Orthogonal dependent variable for 2d fits.

        Returns
        -------
        :func:`tuple`
            A tuple containing the b-spline action matrix; the 'lower' parameter,
            a list of pixel positions, each corresponding to the first
            occurence of position greater than breakpoint indx; and 'upper',
            Same as lower, but denotes the upper pixel positions.
        """
        from .. import uniq
        from ..goddard.math import flegendre
        from .trace import fchebyshev
        nx = x.size
        nbkpt = self.mask.sum()
        if nbkpt < 2*self.nord:
            return (-2, 0, 0)
        n = nbkpt - self.nord
        gb = self.breakpoints[self.mask]
        bw = self.npoly*self.nord
        lower = np.zeros((n - self.nord + 1,), dtype='i4')
        upper = np.zeros((n - self.nord + 1,), dtype='i4') - 1
        indx = self.intrv(x)
        bf1 = self.bsplvn(x, indx)
        action = bf1
        aa = uniq(indx, np.arange(indx.size, dtype='i4'))
        upper[indx[aa]-self.nord+1] = aa
        rindx = indx[::-1]
        bb = uniq(rindx, np.arange(rindx.size, dtype='i4'))
        lower[rindx[bb]-self.nord+1] = nx - bb - 1
        if x2 is not None:
            if x2.size != nx:
                raise ValueError('Dimensions of x and x2 do not match.')
            x2norm = 2.0 * (x2 - self.xmin) / (self.xmax - self.xmin) - 1.0
            if self.funcname == 'poly':
                temppoly = np.ones((nx, self.npoly), dtype='f')
                for i in range(1, self.npoly):
                    temppoly[:, i] = temppoly[:, i-1] * x2norm
            elif self.funcname == 'poly1':
                temppoly = np.tile(x2norm, self.npoly).reshape(nx, self.npoly)
                for i in range(1, self.npoly):
                    temppoly[:, i] = temppoly[:, i-1] * x2norm
            elif self.funcname == 'chebyshev':
                temppoly = fchebyshev(x2norm, self.npoly)
            elif self.funcname == 'legendre':
                temppoly = flegendre(x2norm, self.npoly)
            else:
                raise ValueError('Unknown value of funcname.')
            action = np.zeros((nx, bw), dtype='d')
            counter = -1
            for ii in range(self.nord):
                for jj in range(self.npoly):
                    counter += 1
                    action[:, counter] = bf1[:, ii]*temppoly[:, jj]
        return (action, lower, upper)

    def intrv(self, x):
        """Find the segment between breakpoints which contain each value in the array x.

        The minimum breakpoint is nbkptord -1, and the maximum
        is nbkpt - nbkptord - 1.

        Parameters
        ----------
        x : :class:`numpy.ndarray`
            Data values, assumed to be monotonically increasing.

        Returns
        -------
        :class:`numpy.ndarray`
            Position of array elements with respect to breakpoints.
        """
        gb = self.breakpoints[self.mask]
        n = gb.size - self.nord
        indx = np.zeros((x.size,), dtype='i4')
        ileft = self.nord - 1
        for i in range(x.size):
            while x[i] > gb[ileft+1] and ileft < n - 1:
                ileft += 1
            indx[i] = ileft
        return indx

    def bsplvn(self, x, ileft):
        """To be documented.

        Parameters
        ----------
        x : :class:`numpy.ndarray`
            To be documented.
        ileft : :class:`int`
            To be documented

        Returns
        -------
        :class:`numpy.ndarray`
            To be documented.
        """
        bkpt = self.breakpoints[self.mask]
        vnikx = np.zeros((x.size, self.nord), dtype=x.dtype)
        deltap = vnikx.copy()
        deltam = vnikx.copy()
        j = 0
        vnikx[:, 0] = 1.0
        while j < self.nord - 1:
            ipj = ileft+j+1
            deltap[:, j] = bkpt[ipj] - x
            imj = ileft-j
            deltam[:, j] = x - bkpt[imj]
            vmprev = 0.0
            for l in range(j+1):
                vm = vnikx[:, l]/(deltap[:, l] + deltam[:, j-l])
                vnikx[:, l] = vm*deltap[:, l] + vmprev
                vmprev = vm*deltam[:, j-l]
            j += 1
            vnikx[:, j] = vmprev
        return vnikx

    def value(self, x, x2=None, action=None, lower=None, upper=None):
        """Evaluate a bspline at specified values.

        Parameters
        ----------
        x : :class:`numpy.ndarray`
            Independent variable.
        x2 : :class:`numpy.ndarray`, optional
            Orthogonal dependent variable for 2d fits.
        action : :class:`numpy.ndarray`, optional
            Action matrix to use.  If not supplied it is calculated.
        lower : :class:`numpy.ndarray`, optional
            If the action parameter is supplied, this parameter must also
            be supplied.
        upper : :class:`numpy.ndarray`, optional
            If the action parameter is supplied, this parameter must also
            be supplied.

        Returns
        -------
        :func:`tuple`
            A tuple containing the results of the bspline evaluation and a
            mask indicating where the evaluation was good.
        """
        xsort = x.argsort()
        xwork = x[xsort]
        if x2 is not None:
            x2work = x2[xsort]
        else:
            x2work = None
        if action is not None:
            if lower is None or upper is None:
                raise ValueError('Must specify lower and upper if action is set.')
        else:
            action, lower, upper = self.action(xwork, x2=x2work)
        yfit = np.zeros(x.shape, dtype=x.dtype)
        bw = self.npoly * self.nord
        spot = np.arange(bw, dtype='i4')
        goodbk = self.mask.nonzero()[0]
        coeffbk = self.mask[self.nord:].nonzero()[0]
        n = self.mask.sum() - self.nord
        if self.npoly > 1:
            goodcoeff = self.coeff[:, coeffbk]
        else:
            goodcoeff = self.coeff[coeffbk]
        # maskthis = np.zeros(xwork.shape,dtype=xwork.dtype)
        for i in range(n-self.nord+1):
            ict = upper[i] - lower[i] + 1
            if ict > 0:
                yfit[lower[i]:upper[i]+1] = np.dot(
                    action[lower[i]:upper[i]+1, :], goodcoeff[i*self.npoly+spot])
        yy = yfit.copy()
        yy[xsort] = yfit
        mask = np.ones(x.shape, dtype='bool')
        gb = self.breakpoints[goodbk]
        outside = ((x < gb[self.nord-1]) | (x > gb[n]))
        if outside.any():
            mask[outside] = False
        hmm = ((np.diff(goodbk) > 2).nonzero())[0]
        for jj in range(hmm.size):
            inside = ((x >= self.breakpoints[goodbk[hmm[jj]]]) &
                      (x <= self.breakpoints[goodbk[hmm[jj]+1]-1]))
            if inside.any():
                mask[inside] = False
        return (yy, mask)

    def maskpoints(self, err):
        """Perform simple logic of which breakpoints to mask.

        Parameters
        ----------
        err : :class:`numpy.ndarray`
            The list of indexes returned by the cholesky routines.

        Returns
        -------
        :class:`int`
            An integer indicating the results of the masking.  -1 indicates
            that the error points were successfully masked.  -2 indicates
            failure; the calculation should be aborted.

        Notes
        -----
        The mask attribute is modified, assuming it is possible to create the
        mask.
        """
        nbkpt = self.mask.sum()
        if nbkpt <= 2*self.nord:
            return -2
        hmm = err[np.unique(err/self.npoly)]/self.npoly
        n = nbkpt - self.nord
        if np.any(hmm >= n):
            return -2
        test = np.zeros(nbkpt, dtype='bool')
        for jj in range(-np.ceil(nord/2.0), nord/2.0):
            foo = np.where((hmm+jj) > 0, hmm+jj, np.zeros(hmm.shape, dtype=hmm.dtype))
            inside = np.where((foo+nord) < n-1, foo+nord, np.zeros(hmm.shape, dtype=hmm.dtype)+n-1)
            test[inside] = True
        if test.any():
            reality = self.mask[test]
            if self.mask[reality].any():
                self.mask[reality] = False
                return -1
            else:
                return -2
        else:
            return -2


def cholesky_band(l, mininf=0.0, verbose=False):
    """Compute Cholesky decomposition of banded matrix.

    Parameters
    ----------
    l : :class:`numpy.ndarray`
        A matrix on which to perform the Cholesky decomposition.
    mininf : :class:`float`, optional
        Entries in the `l` matrix are considered negative if they are less
        than this value (default 0.0).
    verbose : :class:`bool`, optional
        If set to ``True``, print some debugging information.

    Returns
    -------
    :func:`tuple`
        If problems were detected, the first item will be the index or
        indexes where the problem was detected, and the second item will simply
        be the input matrix.  If no problems were detected, the first item
        will be -1, and the second item will be the Cholesky decomposition.
    """
    from warnings import warn
    from . import PydlutilsUserWarning
    lower = l.copy()
    bw, nn = lower.shape
    n = nn - bw
    negative = lower[0, 0:n] <= mininf
    if negative.any() or not np.all(np.isfinite(lower)):
        warn('Bad entries: ' + str(negative.nonzero()[0]), PydlutilsUserWarning)
        return (negative.nonzero()[0], l)
    kn = bw - 1
    spot = np.arange(kn, dtype='i4') + 1
    bi = np.arange(kn, dtype='i4')
    for i in range(1, kn):
        bi = np.append(bi, np.arange(kn-i, dtype='i4') + (kn+1)*i)
    for j in range(n):
        lower[0, j] = np.sqrt(lower[0, j])
        lower[spot, j] /= lower[0, j]
        x = lower[spot, j]
        if not np.all(np.isfinite(x)):
            warn('NaN found in cholesky_band.', PydlutilsUserWarning)
            return (j, l)
        hmm = np.outer(x, x)
        here = bi+(j+1)*bw
        lower.T.flat[here] -= hmm.flat[bi]
    return (-1, lower)


def cholesky_solve(a, bb):
    """Solve the equation Ax=b where A is a Cholesky-banded matrix.

    Parameters
    ----------
    a : :class:`numpy.ndarray`
        :math:`A` in :math:`A x = b`.
    bb : :class:`numpy.ndarray`
        :math:`b` in :math:`A x = b`.

    Returns
    -------
    :func:`tuple`
        A tuple containing the status and the result of the solution.  The
        status is always -1.
    """
    b = bb.copy()
    bw = a.shape[0]
    n = b.shape[0] - bw
    kn = bw - 1
    spot = np.arange(kn, dtype='i4') + 1
    for j in range(n):
        b[j] /= a[0, j]
        b[j+spot] -= b[j]*a[spot, j]
    spot = kn - np.arange(kn, dtype='i4')
    for j in range(n-1, -1, -1):
        b[j] = (b[j] - np.sum(a[spot, j] * b[j+spot]))/a[0, j]
    return (-1, b)


def iterfit(xdata, ydata, invvar=None, upper=5, lower=5, x2=None,
            maxiter=10, **kwargs):
    """Iteratively fit a b-spline set to data, with rejection.

    Parameters
    ----------
    xdata : :class:`numpy.ndarray`
        Independent variable.
    ydata : :class:`numpy.ndarray`
        Dependent variable.
    invvar : :class:`numpy.ndarray`
        Inverse variance of `ydata`.  If not set, it will be calculated based
        on the standard deviation.
    upper : :class:`int` or :class:`float`
        Upper rejection threshold in units of sigma, defaults to 5 sigma.
    lower : :class:`int` or :class:`float`
        Lower rejection threshold in units of sigma, defaults to 5 sigma.
    x2 : :class:`numpy.ndarray`, optional
        Orthogonal dependent variable for 2d fits.
    maxiter : :class:`int`, optional
        Maximum number of rejection iterations, default 10.  Set this to
        zero to disable rejection.

    Returns
    -------
    :func:`tuple`
        A tuple containing the fitted bspline object and an output mask.
    """
    from .math import djs_reject
    nx = xdata.size
    if ydata.size != nx:
        raise ValueError('Dimensions of xdata and ydata do not agree.')
    if invvar is not None:
        if invvar.size != nx:
            raise ValueError('Dimensions of xdata and invvar do not agree.')
    else:
        #
        # This correction to the variance makes it the same
        # as IDL's variance()
        #
        var = ydata.var()*(float(nx)/float(nx-1))
        if var == 0:
            var = 1.0
        invvar = np.ones(ydata.shape, dtype=ydata.dtype)/var
    if x2 is not None:
        if x2.size != nx:
            raise ValueError('Dimensions of xdata and x2 do not agree.')
    yfit = np.zeros(ydata.shape)
    if invvar.size == 1:
        outmask = True
    else:
        outmask = np.ones(invvar.shape, dtype='bool')
    xsort = xdata.argsort()
    maskwork = (outmask & (invvar > 0))[xsort]
    if 'oldset' in kwargs:
        sset = kwargs['oldset']
        sset.mask = True
        sset.coeff = 0
    else:
        if not maskwork.any():
            raise ValueError('No valid data points.')
            # return (None,None)
        if 'fullbkpt' in kwargs:
            fullbkpt = kwargs['fullbkpt']
        else:
            sset = bspline(xdata[xsort[maskwork]], **kwargs)
            if maskwork.sum() < sset.nord:
                print('Number of good data points fewer than nord.')
                return (sset, outmask)
            if x2 is not None:
                if 'xmin' in kwargs:
                    xmin = kwargs['xmin']
                else:
                    xmin = x2.min()
                if 'xmax' in kwargs:
                    xmax = kwargs['xmax']
                else:
                    xmax = x2.max()
                if xmin == xmax:
                    xmax = xmin + 1
                sset.xmin = xmin
                sset.xmax = xmax
                if 'funcname' in kwargs:
                    sset.funcname = kwargs['funcname']
    xwork = xdata[xsort]
    ywork = ydata[xsort]
    invwork = invvar[xsort]
    if x2 is not None:
        x2work = x2[xsort]
    else:
        x2work = None
    iiter = 0
    error = 0
    qdone = -1
    while (error != 0 or qdone == -1) and iiter <= maxiter:
        goodbk = sset.mask.nonzero()[0]
        if maskwork.sum() <= 1 or not sset.mask.any():
            sset.coeff = 0
            iiter = maxiter + 1
        else:
            if 'requiren' in kwargs:
                i = 0
                while xwork[i] < sset.breakpoints[goodbk[sset.nord]] and i < nx-1:
                    i += 1
                ct = 0
                for ileft in range(sset.nord, sset.mask.sum()-sset.nord+1):
                    while (xwork[i] >= sset.breakpoints[goodbk[ileft]] and
                           xwork[i] < sset.breakpoints[goodbk[ileft+1]] and
                           i < nx-1):
                        ct += invwork[i]*maskwork[i] > 0
                        i += 1
                    if ct >= kwargs['requiren']:
                        ct = 0
                    else:
                        sset.mask[goodbk[ileft]] = False
            error, yfit = sset.fit(xwork, ywork, invwork*maskwork,
                                   x2=x2work)
        iiter += 1
        inmask = maskwork
        if error == -2:
            return (sset, outmask)
        elif error == 0:
            maskwork, qdone = djs_reject(ywork, yfit, invvar=invwork,
                                         inmask=inmask, outmask=maskwork,
                                         upper=upper, lower=lower)
        else:
            pass
    outmask[xsort] = maskwork
    temp = yfit
    yfit[xsort] = temp
    return (sset, outmask)
