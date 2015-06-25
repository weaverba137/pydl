# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def xy2traceset(xpos,ypos,**kwargs):
    """Convert from x,y positions to a trace set.

    Parameters
    ----------
    xpos, ypos : array-like
        X,Y positions corresponding as [nx,Ntrace] arrays.
    invvar : array-like, optional
        Inverse variances for fitting.
    func : str, optional
        Function type for fitting; defaults to 'legendre'.
    ncoeff : int, optional
        Number of coefficients to fit.  Defaults to 3.
    xmin, xmax : float, optional
        Explicitly set minimum and maximum values, instead of computing
        them from `xpos`.
    maxiter : int, optional
        Maximum number of rejection iterations; set to 0 for no rejection;
        default to 10.
    inmask : array-like, optional
        Mask set to 1 for good points and 0 for rejected points;
        same dimensions as `xpos`, `ypos`.  Points rejected by `inmask`
        are always rejected from the fits (the rejection is "sticky"),
        and will also be marked as rejected in the outmask attribute.
    ia, inputans, inputfunc : array-like, optional
        These arguments will be passed to :func:`func_fit`.
    xjumplo : float, optional
        x position locating start of an x discontinuity
    xjumphi : float, optional
        x position locating end of that x discontinuity
    xjumpval : float, optional
        magnitude of the discontinuity "jump" between those bounds
        (previous 3 keywords motivated by BOSS 2-phase readout)

    Returns
    -------
    xy2traceset : TraceSet
        A :class:`TraceSet` object.
    """
    from . import TraceSet
    return TraceSet(xpos,ypos,**kwargs)
