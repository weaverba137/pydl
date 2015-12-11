# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""This module corresponds to the trace directory in idlutils.
"""
import numpy as np
from . import PydlutilsException
from ..goddard.math import flegendre


def fchebyshev(x, m):
    """Compute the first `m` Chebyshev polynomials.

    Parameters
    ----------
    x : array-like
        Compute the Chebyshev polynomials at these abscissa values.
    m : :class:`int`
        The number of Chebyshev polynomials to compute.  For example, if
        :math:`m = 3`, :math:`T_0 (x)`, :math:`T_1 (x)` and
        :math:`T_2 (x)` will be computed.

    Returns
    -------
    :class:`numpy.ndarray`
    """
    from scipy.special import chebyt
    if isinstance(x, np.ndarray):
        n = x.size
    else:
        n = 1
    if m < 1:
        raise ValueError('Order of Chebyshev polynomial must be at least 1.')
    try:
        dt = x.dtype
    except AttributeError:
        dt = np.float64
    leg = np.ones((m, n), dtype=dt)
    if m >= 2:
        leg[1, :] = x
    if m >= 3:
        for k in range(2, m):
            leg[k, :] = np.polyval(chebyt(k), x)
    return leg


def fchebyshev_split(x, m):
    """Compute the first `m` Chebyshev polynomials, but modified to allow a
    split in the baseline at :math:`x=0`.  The intent is to allow a model fit
    where a constant term is different for positive and negative `x`.

    Parameters
    ----------
    x : array-like
        Compute the Chebyshev polynomials at these abscissa values.
    m : :class:`int`
        The number of Chebyshev polynomials to compute.  For example, if
        :math:`m = 3`, :math:`T_0 (x)`, :math:`T_1 (x)` and
        :math:`T_2 (x)` will be computed.

    Returns
    -------
    :class:`numpy.ndarray`
    """
    import numpy as np
    if isinstance(x, np.ndarray):
        n = x.size
    else:
        n = 1
    if m < 2:
        raise ValueError('Order of polynomial must be at least 2.')
    try:
        dt = x.dtype
    except AttributeError:
        dt = np.float64
    leg = np.ones((m, n), dtype=dt)
    try:
        leg[0, :] = (x >= 0).astype(x.dtype)
    except AttributeError:
        leg[0, :] = np.double(x >= 0)
    if m > 2:
        leg[2, :] = x
    if m > 3:
        for k in range(3, m):
            leg[k, :] = 2.0 * x * leg[k-1, :] - leg[k-2, :]
    return leg


def fpoly(x, m):
    """Compute the first `m` simple polynomials.

    Parameters
    ----------
    x : array-like
        Compute the simple polynomials at these abscissa values.
    m : :class:`int`
        The number of simple polynomials to compute.  For example, if
        :math:`m = 3`, :math:`x^0`, :math:`x^1` and
        :math:`x^2` will be computed.

    Returns
    -------
    :class:`numpy.ndarray`
    """
    if isinstance(x, np.ndarray):
        n = x.size
    else:
        n = 1
    if m < 1:
        raise ValueError('Order of polynomial must be at least 1.')
    try:
        dt = x.dtype
    except AttributeError:
        dt = np.float64
    leg = np.ones((m, n), dtype=dt)
    if m >= 2:
        leg[1, :] = x
    if m >= 3:
        for k in range(2, m):
            leg[k, :] = leg[k-1, :] * x
    return leg


def func_fit(x, y, ncoeff, invvar=None, function_name='legendre', ia=None,
            inputans=None, inputfunc=None):
    """Fit `x`, `y` positions to a functional form.

    Parameters
    ----------
    x : array-like
        X values (independent variable).
    y : array-like
        Y values (dependent variable).
    ncoeff : :class:`int`
        Number of coefficients to fit.
    invvar : array-like, optional
        Weight values; inverse variance.
    function_name : :class:`str`, optional
        Function name, default 'legendre'.
    ia : array-like, optional
        An array of bool of length `ncoeff` specifying free (``True``) and
        fixed (``False``) parameters.
    inputans : array-like, optional
        An array of values of length `ncoeff` specifying the values of
        the fixed parameters.
    inputfunc : array-like, optional
        Multiply the function fit by these values.

    Returns
    -------
    :func:`tuple` of array-like
        Fit coefficients, length `ncoeff`; fitted values.

    Raises
    ------
    KeyError
        If an invalid function type is selected.
    ValueError
        If input dimensions do not agree.
    """
    if x.shape != y.shape:
        raise ValueError('Dimensions of X and Y do not agree!')
    if invvar is None:
        invvar = np.ones(x.shape, dtype=x.dtype)
    else:
        if invvar.shape != x.shape:
            raise ValueError('Dimensions of X and invvar do not agree!')
    if ia is None:
        ia = np.ones((ncoeff,), dtype=np.bool)
    if not ia.all():
        if inputans is None:
            inputans = np.zeros((ncoeff,), dtype=x.dtype)
    #
    # Select unmasked points
    #
    igood = (invvar > 0).nonzero()[0]
    ngood = len(igood)
    res = np.zeros((ncoeff,), dtype=x.dtype)
    yfit = np.zeros(x.shape, dtype=x.dtype)
    if ngood == 0:
        pass
    elif ngood == 1:
        res[0] = y[igood[0]]
        yfit += y[igood[0]]
    else:
        ncfit = min(ngood, ncoeff)
        function_map = {
            'legendre': flegendre,
            'flegendre': flegendre,
            'chebyshev': fchebyshev,
            'fchebyshev': fchebyshev,
            'chebyshev_split': fchebyshev_split,
            'fchebyshev_split': fchebyshev_split,
            'poly': fpoly,
            'fpoly': fpoly
            }
        try:
            legarr = function_map[function_name](x, ncfit)
        except KeyError:
            raise KeyError('Unknown function type: {0}'.format(function_name))
        if inputfunc is not None:
            if inputfunc.shape != x.shape:
                raise ValueError('Dimensions of X and inputfunc do not agree!')
            legarr *= np.tile(inputfunc, ncfit).reshape(ncfit, x.shape[0])
        yfix = np.zeros(x.shape, dtype=x.dtype)
        nonfix = ia[0:ncfit].nonzero()[0]
        nparams = len(nonfix)
        fixed = (~ia[0:ncfit]).nonzero()[0]
        if len(fixed) > 0:
            yfix = np.dot(legarr.T, inputans * (1 - ia))
            ysub = y - yfix
            finalarr = legarr[nonfix, :]
        else:
            finalarr = legarr
            ysub = y
        # extra2 = finalarr * np.outer(np.ones((nparams,), dtype=x.dtype),
        #                             (invvar > 0))
        extra2 = finalarr * np.outer(np.ones((nparams,), dtype=x.dtype),
                                    invvar)
        alpha = np.dot(finalarr, extra2.T)
        # assert alpha.dtype == x.dtype
        if nparams > 1:
            # beta = np.dot(ysub * (invvar > 0), finalarr.T)
            beta = np.dot(ysub * invvar, finalarr.T)
            assert beta.dtype == x.dtype
            # uu,ww,vv = np.linalg.svd(alpha, full_matrices=False)
            res[nonfix] = np.linalg.solve(alpha, beta)
        else:
            # res[nonfix] = (ysub * (invvar > 0) * finalarr).sum()/alpha
            res[nonfix] = (ysub * invvar * finalarr).sum()/alpha
        if len(fixed) > 0:
            res[fixed] = inputans[fixed]
        yfit = np.dot(legarr.T, res[0:ncfit])
    return (res, yfit)


class TraceSet(object):
    """Implements the idea of a trace set.

    Attributes
    ----------
    func : :class:`str`
        Name of function type used to fit the trace set.
    xmin : float-like
        Minimum x value.
    xmax : float-like
        Maximum x value.
    coeff : array-like
        Coefficients of the trace set fit.
    nTrace : :class:`int`
        Number of traces in the object.
    ncoeff : :class:`int`
        Number of coefficients of the trace set fit.
    xjumplo : float-like
        Jump value, for BOSS readouts.
    xjumphi : float-like
        Jump value, for BOSS readouts.
    xjumpval : float-like
        Jump value, for BOSS readouts.
    outmask : array-like
        When initialized with x,y positions, this contains the rejected
        points.
    yfit : array-like
        When initialized with x,y positions, this contains the fitted y
        values.
    """
    _func_map = {'poly': fpoly, 'legendre': flegendre,
                    'chebyshev': fchebyshev}

    def __init__(self, *args, **kwargs):
        """This class can be initialized either with a set of xy positions,
        or with a trace set HDU from a FITS file.
        """
        from astropy.io.fits.fitsrec import FITS_rec
        from .math import djs_reject
        if len(args) == 1 and isinstance(args[0], FITS_rec):
            #
            # Initialize with FITS data
            #
            self.func = args[0]['FUNC'][0]
            self.xmin = args[0]['XMIN'][0]
            self.xmax = args[0]['XMAX'][0]
            self.coeff = args[0]['COEFF'][0]
            self.nTrace = self.coeff.shape[0]
            self.ncoeff = self.coeff.shape[1]
            if 'XJUMPLO' in args[0].dtype.names:
                self.xjumplo = args[0]['XJUMPLO'][0]
                self.xjumphi = args[0]['XJUMPHI'][0]
                self.xjumpval = args[0]['XJUMPVAL'][0]
            else:
                self.xjumplo = None
                self.xjumphi = None
                self.xjumpval = None
            self.outmask = None
            self.yfit = None
        elif len(args) == 2:
            #
            # Initialize with x, y positions.
            #
            xpos = args[0]
            ypos = args[1]
            self.nTrace = xpos.shape[0]
            if 'invvar' in kwargs:
                invvar = kwargs['invvar']
            else:
                invvar = np.ones(xpos.shape, dtype=xpos.dtype)
            if 'func' in kwargs:
                self.func = kwargs['func']
            else:
                self.func = 'legendre'
            if 'ncoeff' in kwargs:
                self.ncoeff = int(kwargs['ncoeff'])
            else:
                self.ncoeff = 3
            if 'xmin' in kwargs:
                self.xmin = np.float64(kwargs['xmin'])
            else:
                self.xmin = xpos.min()
            if 'xmax' in kwargs:
                self.xmax = np.float64(kwargs['xmax'])
            else:
                self.xmax = xpos.max()
            if 'maxiter' in kwargs:
                maxiter = int(kwargs['maxiter'])
            else:
                maxiter = 10
            if 'inmask' in kwargs:
                inmask = kwargs['inmask']
            else:
                inmask = np.ones(xpos.shape, dtype=np.bool)
            do_jump = False
            if 'xjumplo' in kwargs:
                do_jump = True
                self.xjumplo = np.float64(kwargs['xjumplo'])
            else:
                self.xjumplo = None
            if 'xjumphi' in kwargs:
                self.xjumphi = np.float64(kwargs['xjumphi'])
            else:
                self.xjumphi = None
            if 'xjumpval' in kwargs:
                self.xjumpval = np.float64(kwargs['xjumpval'])
            else:
                self.xjumpval = None
            self.coeff = np.zeros((self.nTrace, self.ncoeff), dtype=xpos.dtype)
            self.outmask = np.zeros(xpos.shape, dtype=np.bool)
            self.yfit = np.zeros(xpos.shape, dtype=xpos.dtype)
            for iTrace in range(self.nTrace):
                xvec = self.xnorm(xpos[iTrace, :], do_jump)
                iIter = 0
                qdone = False
                tempivar = (invvar[iTrace, :] *
                            inmask[iTrace, :].astype(invvar.dtype))
                thismask = tempivar > 0
                while (not qdone) and (iIter <= maxiter):
                    res, ycurfit = func_fit(xvec, ypos[iTrace, :], self.ncoeff,
                        invvar=tempivar, function_name=self.func)
                    thismask, qdone = djs_reject(ypos[iTrace, :], ycurfit,
                                                invvar=tempivar)
                    iIter += 1
                self.yfit[iTrace, :] = ycurfit
                self.coeff[iTrace, :] = res
                self.outmask[iTrace, :] = thismask
        else:
            raise PydlutilsException("Wrong number of arguments to TraceSet!")

    def xy(self, xpos=None, ignore_jump=False):
        """Convert from a trace set to an array of x,y positions.

        Parameters
        ----------
        xpos : array-like, optional
            If provided, evaluate the trace set at these positions.  Otherwise
            the positions will be constructed from the trace set object iself.
        ignore_jump : :class:`bool`, optional
            If ``True``, ignore any jump information in the `tset` object

        Returns
        -------
        :func:`tuple` of array-like
            The x, y positions.
        """
        from .misc import djs_laxisgen
        do_jump = self.has_jump and (not ignore_jump)
        if xpos is None:
            xpos = djs_laxisgen([self.nTrace, self.nx], iaxis=1) + self.xmin
        ypos = np.zeros(xpos.shape, dtype=xpos.dtype)
        for iTrace in range(self.nTrace):
            xvec = self.xnorm(xpos[iTrace, :], do_jump)
            legarr = self._func_map[self.func](xvec, self.ncoeff)
            ypos[iTrace, :] = np.dot(legarr.T, self.coeff[iTrace, :])
        return (xpos, ypos)

    @property
    def has_jump(self):
        """``True`` if jump conditions are set.
        """
        return self.xjumplo is not None

    @property
    def xRange(self):
        """Range of x values.
        """
        return self.xmax - self.xmin

    @property
    def nx(self):
        """Number of x values.
        """
        return int(self.xRange + 1)

    @property
    def xmid(self):
        """Midpoint of x values.
        """
        return 0.5 * (self.xmin + self.xmax)

    def xnorm(self, xinput, jump):
        """Convert input x coordinates to normalized coordinates suitable
        for input to special polynomials.

        Parameters
        ----------
        xinput : array-like
            Input coordinates.
        jump : :class:`bool`
            Set to ``True`` if there is a jump.

        Returns
        -------
        array-like
            Normalized coordinates.
        """
        if jump:
            # Vector specifying what fraction of the jump has passed:
            jfrac = np.minimum(np.maximum(((xinput - self.xjumplo) /
                                (self.xjumphi - self.xjumplo)), 0.), 1.)
            # Conversion to "natural" x baseline:
            xnatural = xinput + jfrac * self.xjumpval
        else:
            xnatural = xinput
        return 2.0 * (xnatural - self.xmid)/self.xRange


def traceset2xy(tset, xpos=None, ignore_jump=False):
    """Convert from a trace set to an array of x,y positions.

    Parameters
    ----------
    tset : :class:`TraceSet`
        A :class:`TraceSet` object.
    xpos : array-like, optional
        If provided, evaluate the trace set at these positions.  Otherwise
        the positions will be constructed from the trace set object iself.
    ignore_jump : bool, optional
        If ``True``, ignore any jump information in the `tset` object

    Returns
    -------
    :func:`tuple` of array-like
        The x, y positions.
    """
    return tset.xy(xpos, ignore_jump)


def xy2traceset(xpos, ypos, **kwargs):
    """Convert from x,y positions to a trace set.

    Parameters
    ----------
    xpos, ypos : array-like
        X,Y positions corresponding as [nx,Ntrace] arrays.
    invvar : array-like, optional
        Inverse variances for fitting.
    func : :class:`str`, optional
        Function type for fitting; defaults to 'legendre'.
    ncoeff : :class:`int`, optional
        Number of coefficients to fit.  Defaults to 3.
    xmin, xmax : :class:`float`, optional
        Explicitly set minimum and maximum values, instead of computing
        them from `xpos`.
    maxiter : :class:`int`, optional
        Maximum number of rejection iterations; set to 0 for no rejection;
        default to 10.
    inmask : array-like, optional
        Mask set to 1 for good points and 0 for rejected points;
        same dimensions as `xpos`, `ypos`.  Points rejected by `inmask`
        are always rejected from the fits (the rejection is "sticky"),
        and will also be marked as rejected in the outmask attribute.
    ia, inputans, inputfunc : array-like, optional
        These arguments will be passed to :func:`func_fit`.
    xjumplo : :class:`float`, optional
        x position locating start of an x discontinuity
    xjumphi : :class:`float`, optional
        x position locating end of that x discontinuity
    xjumpval : :class:`float`, optional
        magnitude of the discontinuity "jump" between those bounds
        (previous 3 keywords motivated by BOSS 2-phase readout)

    Returns
    -------
    :class:`TraceSet`
        A :class:`TraceSet` object.
    """
    return TraceSet(xpos, ypos, **kwargs)
