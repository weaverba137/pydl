#
# $Id$
#
def spherematch(ra1,dec1,ra2,dec2,matchlength,**kwargs):
    """Match points on a sphere.
    """
    import numpy as np
    from pydlutils.spheregroup import chunks, separation
    #
    # Set default values
    #
    if 'chunksize' in kwargs:
        chunksize = kwargs['chunksize']
    else:
        chunksize = max(4.0*matchlength,0.1)
    if 'maxmatch' in kwargs:
        maxmatch = kwargs['maxmatch']
    else:
        maxmatch = 1
    #
    # Initialize chunks
    #
    chunk = chunks(ra1,dec1,chunksize)
    chunk.assign(ra2,dec2,matchlength)
    if 'debug' in kwargs:
        print "raOffset = %7.3f." % chunk.raOffset
        print chunk.raBounds
        print chunk.decBounds
        print chunk.chunkList
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
                sep = separation(ra1[i],dec1[i],ra2[k],dec2[k])
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
