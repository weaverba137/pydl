# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
from astropy.io.fits.fitsrec import FITS_rec
from . import fpoly, fchebyshev, func_fit
from .. import PydlutilsException
from ..math import djs_reject
from ..misc import djs_laxisgen
from ...goddard.math import flegendre
#
class TraceSet(object):
    """Implements the idea of a trace set.

    Attributes
    ----------
    func : str
        Name of function type used to fit the trace set.
    xmin : float-like
        Minimum x value.
    xmax : float-like
        Maximum x value.
    coeff : array-like
        Coefficients of the trace set fit.
    nTrace : int
        Number of traces in the object.
    ncoeff : int
        Number of coefficients of the trace set fit.
    xjumplo : float-like
        Jump value, for BOSS readouts.
    xjumphi : float-like
        Jump value, for BOSS readouts.
    xjumpval : float-like
        Jump value, for BOSS readouts.
    outmask : array-like.
        When initialized with x,y positions, this contains the rejected points.
    yfit : array-like.
        When initialized with x,y positions, this contains the fitted y values.
    """
    _func_map = {'poly':fpoly,'legendre':flegendre,'chebyshev':fchebyshev}
    #
    #
    #
    def __init__(self,*args,**kwargs):
        """This class can be initialized either with a set of xy positions,
        or with a trace set HDU from a FITS file.
        """
        if len(args) == 1 and isinstance(args[0],FITS_rec):
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
                invvar = np.ones(xpos.shape,dtype=xpos.dtype)
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
                inmask = np.ones(xpos.shape,dtype=np.bool)
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
            self.coeff = np.zeros((self.nTrace,self.ncoeff),dtype=xpos.dtype)
            self.outmask = np.zeros(xpos.shape,dtype=np.bool)
            self.yfit = np.zeros(xpos.shape,dtype=xpos.dtype)
            for iTrace in range(self.nTrace):
                xvec = self.xnorm(xpos[iTrace,:],do_jump)
                iIter = 0
                qdone = False
                tempivar = invvar[iTrace,:] * inmask[iTrace,:].astype(invvar.dtype)
                thismask = tempivar > 0
                while (not qdone) and (iIter <= maxiter):
                    res,ycurfit = func_fit(xvec, ypos[iTrace,:], self.ncoeff,
                        invvar=tempivar, function_name=self.func)
                    thismask,qdone = djs_reject(ypos[iTrace,:],ycurfit,invvar=tempivar)
                    iIter += 1
                self.yfit[iTrace,:] = ycurfit
                self.coeff[iTrace,:] = res
                self.outmask[iTrace,:] = thismask
        else:
            raise PydlutilsException("Wrong number of arguments to TraceSet!")
    #
    #
    #
    def xy(self,xpos=None,ignore_jump=False):
        """Convert from a trace set to an array of x,y positions.

        Parameters
        ----------
        xpos : array-like, optional
            If provided, evaluate the trace set at these positions.  Otherwise
            the positions will be constructed from the trace set object iself.
        ignore_jump : bool, optional
            If ``True``, ignore any jump information in the `tset` object

        Returns
        -------
        xy : tuple of array-like
            The x, y positions.
        """
        do_jump = self.has_jump and (not ignore_jump)
        if xpos is None:
            xpos = djs_laxisgen([self.nTrace,self.nx], iaxis=1) + self.xmin
        ypos = np.zeros(xpos.shape,dtype=xpos.dtype)
        for iTrace in range(self.nTrace):
            xvec = self.xnorm(xpos[iTrace,:],do_jump)
            legarr = self._func_map[self.func](xvec,self.ncoeff)
            ypos[iTrace:] = np.dot(legarr.T,self.coeff[iTrace,:])
        return (xpos,ypos)
    #
    #
    #
    @property
    def has_jump(self):
        return self.xjumplo is not None
    #
    #
    #
    @property
    def xRange(self):
        return self.xmax - self.xmin
    #
    #
    #
    @property
    def nx(self):
        return int(self.xRange + 1)
    #
    #
    #
    @property
    def xmid(self):
        return 0.5 * (self.xmin + self.xmax)
    #
    #
    #
    def xnorm(self,xinput,jump):
        """Convert input x coordinates to normalized coordinates suitable
        for input to special polynomials.
        """
        if jump:
            # Vector specifying what fraction of the jump has passed:
            jfrac = np.minimum(np.maximum(((xinput - self.xjumplo) / (self.xjumphi - self.xjumplo)),0.),1.)
            # Conversion to "natural" x baseline:
            xnatural = xinput + jfrac * self.xjumpval
        else:
            xnatural = xinput
        return 2.0 * (xnatural-self.xmid)/self.xRange
