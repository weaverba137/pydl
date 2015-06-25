# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
from . import fpoly, fchebyshev
from ..misc import djs_laxisgen
from ...goddard.math import flegendre
#
class TraceSet(object):
    """Implements the idea of a trace set.

    Attributes
    ----------
    func : str
        Name of function type used to fit the trace set.
    ncoeff : int
        Number of coefficients to fit.
    """
    _func_map = {'poly':fpoly,'legendre':flegendre,'chebyshev':fchebyshev}
    #
    #
    #
    def __init__(self,xpos,ypos,invvar=None,func='legendre',ncoeff=3,
        xmin=None,xmax=None,maxiter=10,inmask=None,
        ia=None,inputans=None,inputfunc=None,
        xjumplo=None,xjumphi=None,xjumpval=None):
        #
        # Capture inputs.
        #
        self.func = func
        self.ncoeff = ncoeff
        if xmin is None:
            self.xmin = xpos.min()
        else:
            self.xmin = xmin
        if xmax is None:
            self.xmax = xpos.max()
        else:
            self.xmax = xmax
        self.maxiter = maxiter
        self.inmask = inmask
        self.xjumplo = xjumplo
        self.xjumphi = xjumphi
        self.xjumpval = xjumpval
        self.ndim = len(ypos.shape)
        self.nx = ypos.shape[0]
        try:
            self.nTrace = ypos.shape[1]
        except IndexError:
            self.nTrace = 1
        self.outmask = np.zeros((self.nx,self.nTrace),dtype=np.bool)
        return
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
        do_jump = (self.xjumpval is not None) and (not ignore_jump)
        xRange = self.xmax - self.xmin
        nx = int(xRange + 1)
        xmid = 0.5 * (self.xmin + self.xmax)
        if xpos is None:
            xpos = djs_laxisgen([nx, self.nTrace], iaxis=0) + self.xmin
        ypos = np.zeros(xpos.shape,dtype=xpos.dtype)
        for iTrace in range(self.nTrace):
            xinput = xpos[:,iTrace]
            if do_jump:
                # Vector specifying what fraction of the jump has passed:
                jfrac = np.minimum(np.maximum(((xinput - self.xjumplo) / (self.xjumphi - self.xjumplo)),0.),1.)
                # Conversion to "natural" x baseline:
                xnatural = xinput + jfrac * self.xjumpval
            else:
                xnatural = xinput
            xvec = 2.0 * (xnatural-xmid)/xRange
            legarr = self._func_map[self.func](xvec,ncoeff)
            ypos[:,iTrace] = np.dot(legarr,self.coeff[:,iTrace])
        return (xpos,ypos)
