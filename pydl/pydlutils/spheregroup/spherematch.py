# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from __future__ import print_function
#
def spherematch(ra1,dec1,ra2,dec2,matchlength,chunksize=None,maxmatch=1,debug=False):
    """Match points on a sphere.

    Parameters
    ----------
    ra1, dec1, ra2, dec2 : ndarray
        The sets of coordinates to match.  Assumed to be in decimal degrees
    matchlength : float
        Two points closer than this separation are matched. Assumed to be in decimal degrees.
    chunksize : float, optional
        Value to pass to chunk assignment.
    maxmatch : int, optional
        Allow up to `maxmatch` matches per coordinate.  Default 1. If set to zero,
        All possible matches will be returned.
    debug : bool, optional

    Returns
    -------
    """
    import numpy as np
    from . import chunks
    from ...goddard.astro import gcirc
    #
    # Set default values
    #
    if chunksize is None:
        chunksize = max(4.0*matchlength,0.1)
    #
    # Initialize chunks
    #
    chunk = chunks(ra1,dec1,chunksize)
    chunk.assign(ra2,dec2,matchlength)
    if debug:
        print("raOffset = {0:7.3f}.".format(chunk.raOffset))
        print(chunk.raBounds)
        print(chunk.decBounds)
        print(chunk.chunkList)
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
