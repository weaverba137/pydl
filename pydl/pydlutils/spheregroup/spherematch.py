# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from __future__ import print_function
#
def spherematch(ra1,dec1,ra2,dec2,matchlength,chunksize=None,maxmatch=1):
    """Match points on a sphere.

    Parameters
    ----------
    ra1, dec1, ra2, dec2 : :class:`numpy.ndarray`
        The sets of coordinates to match.  Assumed to be in decimal degrees
    matchlength : :class:`float`
        Two points closer than this separation are matched. Assumed to be in decimal degrees.
    chunksize : :class:`float`, optional
        Value to pass to chunk assignment.
    maxmatch : :class:`int`, optional
        Allow up to `maxmatch` matches per coordinate.  Default 1. If set to zero,
        All possible matches will be returned.

    Returns
    -------
    spherematch : :func:`tuple`
        A tuple containing the indices into the first set of points, the
        indices into the second set of points and the match distance in
        decimal degrees.

    Notes
    -----
    If you have sets of coordinates that differ in size, call this function
    with the larger list first.  This exploits the inherent asymmetry in the
    underlying code to reduce memory use.

    .. warning:: Behavior at the poles is not well tested.
    """
    import numpy as np
    from . import chunks
    from .. import PydlutilsException
    from ...goddard.astro import gcirc
    #
    # Set default values
    #
    if chunksize is None:
        chunksize = max(4.0*matchlength,0.1)
    #
    # Check input size
    #
    if ra1.size == 1:
        raise PydlutilsException("Change the order of the sets of coordinates!")
    #
    # Initialize chunks
    #
    chunk = chunks(ra1,dec1,chunksize)
    chunk.assign(ra2,dec2,matchlength)
    #
    # Create return arrays
    #
    match1 = list()
    match2 = list()
    distance12 = list()
    for i in range(ra1.size):
        currra = np.fmod(ra1[i]+chunk.raOffset,360.0)
        rachunk,decchunk = chunk.get(currra,dec1[i])
        jmax = len(chunk.chunkList[decchunk][rachunk])
        if jmax > 0:
            for j in range(jmax):
                k = chunk.chunkList[decchunk][rachunk][j]
                sep = gcirc(ra1[i],dec1[i],ra2[k],dec2[k],units=2)/3600.0
                if sep < matchlength:
                    match1.append(i)
                    match2.append(k)
                    distance12.append(sep)
    #
    # Sort distances
    #
    omatch1 = np.array(match1)
    omatch2 = np.array(match2)
    odistance12 = np.array(distance12)
    s = odistance12.argsort()
    #
    # Retain only desired matches
    #
    if maxmatch > 0:
        gotten1 = np.zeros(ra1.size,dtype='i4')
        gotten2 = np.zeros(ra2.size,dtype='i4')
        nmatch = 0
        for i in range(omatch1.size):
            if (gotten1[omatch1[s[i]]] < maxmatch and
                gotten2[omatch2[s[i]]] < maxmatch):
                gotten1[omatch1[s[i]]] += 1
                gotten2[omatch2[s[i]]] += 1
                nmatch += 1
        match1 = np.zeros(nmatch,dtype='i4')
        match2 = np.zeros(nmatch,dtype='i4')
        distance12 = np.zeros(nmatch,dtype='d')
        gotten1[:] = 0
        gotten2[:] = 0
        nmatch = 0
        for i in range(omatch1.size):
            if (gotten1[omatch1[s[i]]] < maxmatch and
                gotten2[omatch2[s[i]]] < maxmatch):
                gotten1[omatch1[s[i]]] += 1
                gotten2[omatch2[s[i]]] += 1
                match1[nmatch] = omatch1[s[i]]
                match2[nmatch] = omatch2[s[i]]
                distance12[nmatch] = odistance12[s[i]]
                nmatch += 1
    else:
        match1 = omatch1[s]
        match2 = omatch2[s]
        distance12 = odistance12[s]
    return (match1,match2,distance12)
