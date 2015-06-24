# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def traceset2xy(tset,xpos=None,ignore_jump=False):
    """Convert from a trace set to an array of x,y positions.

    Parameters
    ----------
    tset : TraceSet
        A :class:`TraceSet` object.
    xpos : array-like, optional
        If provided, evaluate the trace set at these positions.  Otherwise
        the positions will be constructed from the trace set object iself.
    ignore_jump : bool, optional
        If ``True``, ignore any jump information in the `tset` object

    Returns
    -------
    traceset2xy : tuple of array-like
        The x, y positions.
    """
    return tset.xy(xpos,ignore_jump)
