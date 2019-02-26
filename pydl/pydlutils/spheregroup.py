# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""This module corresponds to the spheregroup directory in idlutils.
"""
from warnings import warn
import numpy as np
from . import PydlutilsException, PydlutilsUserWarning
from ..goddard.astro import gcirc


class chunks(object):
    """chunks class

    Functions for creating and manipulating spherical chunks are implemented
    as methods on this class.
    """

    def __init__(self, ra, dec, minSize):
        """Init creates an object whose attributes are similar those created
        by the setchunks() function in the spheregroup library.
        """
        #
        # Save the value of minSize
        #
        self.minSize = minSize
        #
        # Find maximum and minimum dec (in degrees)
        #
        decMin = dec.min()
        decMax = dec.max()
        decRange = decMax - decMin
        #
        # Find the declination boundaries; make them an integer multiple of
        # minSize, with extra room (one cell) on the edges.
        #
        self.nDec = 3 + int(np.floor(decRange/minSize))
        decRange = minSize*float(self.nDec)
        decMin = decMin - 0.5*(decRange - decMax + decMin)
        decMax = decMin + decRange
        if decMin < -90.0 + 3.0*minSize:
            decMin = -90.0
        if decMax > 90.0 - 3.0*minSize:
            decMax = 90.0
        self.decBounds = decMin + ((decMax - decMin) * np.arange(self.nDec + 1,
                                    dtype='d'))/float(self.nDec)
        #
        # Find ra offset which minimizes the range in ra (this should take care
        # of the case that ra crosses zero in some parts
        #
        if abs(self.decBounds[self.nDec]) > abs(self.decBounds[0]):
            cosDecMin = np.cos(np.deg2rad(self.decBounds[self.nDec]))
        else:
            cosDecMin = np.cos(np.deg2rad(self.decBounds[0]))
        if cosDecMin <= 0.0:
            raise PydlutilsException("cosDecMin={0:f} not positive in setchunks().".format(cosDecMin))
        self.raRange, self.raOffset = self.rarange(ra, minSize/cosDecMin)
        self.raMin, self.raMax = self.getraminmax(ra, self.raOffset)
        #
        # Isn't this redundant?
        #
        self.raRange = self.raMax - self.raMin
        #
        # For each declination slice, find the number of ra divisions
        # necessary and set them
        #
        self.raBounds = list()
        self.nRa = list()
        for i in range(self.nDec):
            #
            # Get maximum declination and its cosine
            #
            if abs(self.decBounds[i]) > abs(self.decBounds[i+1]):
                cosDecMin = np.cos(np.deg2rad(self.decBounds[i]))
            else:
                cosDecMin = np.cos(np.deg2rad(self.decBounds[i+1]))
            if cosDecMin <= 0.0:
                raise PydlutilsException("cosDecMin={0:f} not positive in setchunks().".format(cosDecMin))
            #
            # Get raBounds array for this declination array, leave an extra
            # cell on each end
            #
            self.nRa.append(3 + int(np.floor(cosDecMin*self.raRange/minSize)))
            raRangeTmp = minSize*float(self.nRa[i])/cosDecMin
            raMinTmp = self.raMin - 0.5*(raRangeTmp-self.raMax+self.raMin)
            raMaxTmp = raMinTmp + raRangeTmp
            #
            # If we cannot avoid the 0/360 point, embrace it
            #
            if (raRangeTmp >= 360.0 or
                    raMinTmp <= minSize/cosDecMin or
                    raMaxTmp >= 360.0 - minSize/cosDecMin or
                    abs(self.decBounds[i]) == 90.0):
                raMinTmp = 0.0
                raMaxTmp = 360.0
                raRangeTmp = 360.0
            if self.decBounds[i] == -90.0 or self.decBounds[i+1] == 90.0:
                self.nRa[i] = 1
            self.raBounds.append(raMinTmp +
                (raMaxTmp - raMinTmp) * np.arange(self.nRa[i] + 1, dtype='d') /
                float(self.nRa[i]))
        #
        # Create an empty set of lists to hold the output of self.assign()
        #
        self.chunkList = [[list() for j in range(self.nRa[i])] for i in range(self.nDec)]
        #
        # nChunkMax will be the length of the largest list in chunkList
        # it is computed by chunks.assign()
        #
        self.nChunkMax = 0
        return

    def rarange(self, ra, minSize):
        """Finds the offset which yields the smallest raRange & returns both.

        Notes
        -----

        .. warning:: This is not (yet) well-defined for the case of only one point.
        """
        NRA = 6
        raRangeMin = 361.0
        raOffset = 0.0
        EPS = 1.0e-5
        for j in range(NRA):
            raMin, raMax = self.getraminmax(ra, 360.0*float(j)/float(NRA))
            raRange = raMax-raMin
            if (2.0*(raRange-raRangeMin)/(raRange+raRangeMin) < -EPS and
                    raMin > minSize and raMax < 360.0 - minSize):
                raRangeMin = raRange
                raOffset = 360.0*float(j)/float(NRA)
        return (raRangeMin, raOffset)

    def getraminmax(self, ra, raOffset):
        """Utility function used by rarange.
        """
        currRa = np.fmod(ra + raOffset, 360.0)
        return (currRa.min(), currRa.max())

    def cosDecMin(self, i):
        """Frequently used utility function.
        """
        if abs(self.decBounds[i]) > abs(self.decBounds[i+1]):
            return np.cos(np.deg2rad(self.decBounds[i]))
        else:
            return np.cos(np.deg2rad(self.decBounds[i+1]))

    def assign(self, ra, dec, marginSize):
        """Take the objects and the chunks (already defined in the constructor)
        and assign the objects to the appropriate chunks, with some leeway
        given by the parameter marginSize.  Basically, at the end, each
        chunk should be associated with a list of the objects that belong
        to it.
        """
        if marginSize >= self.minSize:
            raise PydlutilsException("marginSize>=minSize ({0:f}={1:f}) in chunks.assign().".format(marginSize, self.minSize))
        chunkDone = [[False for j in range(self.nRa[i])] for i in range(self.nDec)]
        for i in range(ra.size):
            currRa = np.fmod(ra[i] + self.raOffset, 360.0)
            try:
                raChunkMin, raChunkMax, decChunkMin, decChunkMax = self.getbounds(currRa, dec[i], marginSize)
            except PydlutilsException:
                continue
            #
            # Reset chunkDone.  This is silly, but is necessary to
            # reproduce the logic.
            #
            for decChunk in range(decChunkMin, decChunkMax+1):
                for raChunk in range(raChunkMin[decChunk-decChunkMin]-1, raChunkMax[decChunk-decChunkMin]+2):
                    if raChunk < 0:
                        currRaChunk = (raChunk+self.nRa[decChunk]) % self.nRa[decChunk]
                    elif raChunk > self.nRa[decChunk]-1:
                        currRaChunk = (raChunk-self.nRa[decChunk]) % self.nRa[decChunk]
                    else:
                        currRaChunk = raChunk
                    if currRaChunk >= 0 and currRaChunk <= self.nRa[decChunk]-1:
                        chunkDone[decChunk][currRaChunk] = False
            for decChunk in range(decChunkMin, decChunkMax+1):
                for raChunk in range(raChunkMin[decChunk-decChunkMin], raChunkMax[decChunk-decChunkMin]+1):
                    if raChunk < 0:
                        currRaChunk = (raChunk+self.nRa[decChunk]) % self.nRa[decChunk]
                    elif raChunk > self.nRa[decChunk]-1:
                        currRaChunk = (raChunk-self.nRa[decChunk]) % self.nRa[decChunk]
                    else:
                        currRaChunk = raChunk
                    if currRaChunk >= 0 and currRaChunk <= self.nRa[decChunk]-1:
                        if not chunkDone[decChunk][currRaChunk]:
                            self.chunkList[decChunk][currRaChunk].append(i)
                            #
                            # Update nChunkMax
                            #
                            if len(self.chunkList[decChunk][currRaChunk]) > self.nChunkMax:
                                self.nChunkMax = len(self.chunkList[decChunk][currRaChunk])
                            chunkDone[decChunk][currRaChunk] = True
        return

    def getbounds(self, ra, dec, marginSize):
        """Find the set of chunks a point (with margin) belongs to.
        """
        #
        # Find the declination slice without regard to marginSize
        #
        decChunkMin = int(np.floor((dec - self.decBounds[0]) *
            float(self.nDec) /
            (self.decBounds[self.nDec]-self.decBounds[0])))
        decChunkMax = decChunkMin
        if decChunkMin < 0 or decChunkMin > self.nDec - 1:
            raise PydlutilsException("decChunkMin out of range in chunks.getbounds().")
        #
        # Set minimum and maximum bounds of dec
        #
        while dec - self.decBounds[decChunkMin] < marginSize and decChunkMin > 0:
            decChunkMin -= 1
        while self.decBounds[decChunkMax+1] - dec < marginSize and decChunkMax < self.nDec - 1:
            decChunkMax += 1
        #
        # Find ra chunk bounds for each dec chunk
        #
        raChunkMin = np.zeros(decChunkMax-decChunkMin+1, dtype='i4')
        raChunkMax = np.zeros(decChunkMax-decChunkMin+1, dtype='i4')
        for i in range(decChunkMin, decChunkMax+1):
            cosDecMin = self.cosDecMin(i)
            raChunkMin[i-decChunkMin] = int(np.floor((ra - self.raBounds[i][0]) *
                float(self.nRa[i]) /
                (self.raBounds[i][self.nRa[i]] - self.raBounds[i][0])))
            raChunkMax[i-decChunkMin] = raChunkMin[i-decChunkMin]
            if raChunkMin[i-decChunkMin] < 0 or raChunkMin[i-decChunkMin] > self.nRa[i]-1:
                raise PydlutilsException("raChunkMin out of range in chunks.getbounds().")
            #
            # Set minimum and maximum bounds of ra
            #
            raCheck = raChunkMin[i-decChunkMin]
            keepGoing = True
            while keepGoing and raCheck > -1:
                if raCheck >= 0 and raCheck < self.nRa[i]:
                    keepGoing = (ra - self.raBounds[i][raCheck])*cosDecMin < marginSize
                else:
                    keepGoing = False
                if keepGoing:
                    raCheck -= 1
            raChunkMin[i-decChunkMin] = raCheck
            raCheck = raChunkMax[i-decChunkMin]
            keepGoing = True
            while keepGoing and raCheck < self.nRa[i]:
                if raCheck >= 0 and raCheck < self.nRa[i]:
                    keepGoing = (self.raBounds[i][raCheck+1]-ra)*cosDecMin < marginSize
                else:
                    keepGoing = False
                if keepGoing:
                    raCheck += 1
            raChunkMax[i-decChunkMin] = raCheck
        return (raChunkMin, raChunkMax, decChunkMin, decChunkMax)

    def get(self, ra, dec):
        """Find the chunk to which a given point belongs.
        """
        #
        # Find dec chunk
        #
        decChunk = int(np.floor((dec - self.decBounds[0]) *
            float(self.nDec) /
            (self.decBounds[self.nDec]-self.decBounds[0])))
        #
        # Find ra chunk
        #
        if decChunk < self.nDec and decChunk >= 0:
            raChunk = int(np.floor((ra - self.raBounds[decChunk][0]) *
                float(self.nRa[decChunk]) /
                (self.raBounds[decChunk][self.nRa[decChunk]] - self.raBounds[decChunk][0])))
            if raChunk < 0 or raChunk > self.nRa[decChunk]-1:
                raise PydlutilsException("raChunk out of range in chunks.get()")
        else:
            raChunk = -1
        return (raChunk, decChunk)

    def friendsoffriends(self, ra, dec, linkSep):
        """Friends-of-friends using chunked data.
        """
        nPoints = ra.size
        inGroup = np.zeros(nPoints, dtype='i4') - 1
        #
        # mapGroups contains an equivalency mapping of groups.  mapGroup[i]=j
        # means i and j are actually the same group.  j<=i always, by design.
        # The largest number of groups you can get
        # (assuming linkSep < marginSize < minSize) is 9 times the number of
        # targets
        #
        mapGroups = np.zeros(9*nPoints, dtype='i4') - 1
        nMapGroups = 0
        for i in range(self.nDec):
            for j in range(self.nRa[i]):
                if len(self.chunkList[i][j]) > 0:
                    chunkGroup = self.chunkfriendsoffriends(ra, dec, self.chunkList[i][j], linkSep)
                    for k in range(chunkGroup.nGroups):
                        minEarly = 9*nPoints
                        l = chunkGroup.firstGroup[k]
                        while l != -1:
                            if inGroup[self.chunkList[i][j][l]] != -1:
                                checkEarly = inGroup[self.chunkList[i][j][l]]
                                while mapGroups[checkEarly] != checkEarly:
                                    checkEarly = mapGroups[checkEarly]
                                minEarly = min(minEarly, checkEarly)
                            else:
                                inGroup[self.chunkList[i][j][l]] = nMapGroups
                            l = chunkGroup.nextGroup[l]
                        if minEarly == 9*nPoints:
                            mapGroups[nMapGroups] = nMapGroups
                        else:
                            mapGroups[nMapGroups] = minEarly
                            l = chunkGroup.firstGroup[k]
                            while l != -1:
                                checkEarly = inGroup[self.chunkList[i][j][l]]
                                while mapGroups[checkEarly] != checkEarly:
                                    tmpEarly = mapGroups[checkEarly]
                                    mapGroups[checkEarly] = minEarly
                                    checkEarly = tmpEarly
                                mapGroups[checkEarly] = minEarly
                                l = chunkGroup.nextGroup[l]
                        nMapGroups += 1
        #
        # Now all groups which are mapped to themselves are the real groups
        # Make sure the mappings are set up to go all the way down.
        #
        nGroups = 0
        for i in range(nMapGroups):
            if mapGroups[i] != -1:
                if mapGroups[i] == i:
                    mapGroups[i] = nGroups
                    nGroups += 1
                else:
                    mapGroups[i] = mapGroups[mapGroups[i]]
            else:
                raise PydlutilsException("MapGroups[{0:d}]={1:d} in chunks.friendsoffriends().".format(i, mapGroups[i]))
        for i in range(nPoints):
            inGroup[i] = mapGroups[inGroup[i]]
        firstGroup = np.zeros(nPoints, dtype='i4') - 1
        nextGroup = np.zeros(nPoints, dtype='i4') - 1
        multGroup = np.zeros(nPoints, dtype='i4')
        for i in range(nPoints-1, -1, -1):
            nextGroup[i] = firstGroup[inGroup[i]]
            firstGroup[inGroup[i]] = i
        for i in range(nGroups):
            j = firstGroup[i]
            while j != -1:
                multGroup[i] += 1
                j = nextGroup[j]
        return (inGroup, multGroup, firstGroup, nextGroup, nGroups)

    def chunkfriendsoffriends(self, ra, dec, chunkList, linkSep):
        """Does friends-of-friends on the ra, dec that are defined by
        chunkList.
        """
        #
        # Convert ra, dec into something that can be digested by the
        # groups object.
        #
        x = np.deg2rad(np.vstack((ra[chunkList], dec[chunkList])))
        radLinkSep = np.deg2rad(linkSep)
        group = groups(x, radLinkSep, 'sphereradec')
        return group


class groups(object):
    """Group a set of objects (a list of coordinates in some space) based on
    a friends-of-friends algorithm
    """

    @staticmethod
    def euclid(x1, x2):
        """Pythagorean theorem in Euclidean space with arbitrary number
        of dimensions.
        """
        return np.sqrt(((x1-x2)**2).sum())

    @staticmethod
    def sphereradec(x1, x2):
        """Separation of two points on a 2D-sphere, assuming they are in
        longitude-latitude or right ascension-declination form.  Assumes
        everything is already in radians.
        """
        return gcirc(x1[0], x1[1], x2[0], x2[1], units=0)

    def __init__(self, coordinates, distance, separation='euclid'):
        """Init creates an object and performs the friends-of-friends
        algorithm.  The coordinates can have arbitrary dimensions, with each
        column representing one of the dimensions.  Each row defines an object.
        If separation is not defined it defaults to Euclidean space.
        """
        #
        # Find a separation function
        #
        if callable(separation):
            self.separation = separation
        elif isinstance(separation, (str,)):
            if separation == 'euclid':
                self.separation = self.euclid
            elif separation == 'sphereradec':
                self.separation = self.sphereradec
            else:
                raise PydlutilsException("Unknown separation function: {0}.".format(separation))
        else:
            raise PydlutilsException("Improper type for separation!")
        #
        # Save information about the coordinates.
        #
        nGroups = 0
        nTargets = coordinates.shape[1]
        multGroup = np.zeros(nTargets, dtype='i4')
        firstGroup = np.zeros(nTargets, dtype='i4') - 1
        nextGroup = np.zeros(nTargets, dtype='i4') - 1
        inGroup = np.arange(nTargets, dtype='i4')
        #
        # Find all the other targets associated with each target
        #
        for i in range(nTargets):
            nTmp = 0
            minGroup = nGroups
            for j in range(nTargets):
                sep = self.separation(coordinates[:, i], coordinates[:, j])
                if sep <= distance:
                    multGroup[nTmp] = j
                    minGroup = min(minGroup, inGroup[j])
                    nTmp += 1
            #
            # Use this minimum for all
            #
            for j in range(nTmp):
                if inGroup[multGroup[j]] < nTargets:
                    k = firstGroup[inGroup[multGroup[j]]]
                    while k != -1:
                        inGroup[k] = minGroup
                        k = nextGroup[k]
                inGroup[multGroup[j]] = minGroup
            #
            # If it is a new group (no earlier groups), increment nGroups
            #
            if minGroup == nGroups:
                nGroups += 1
            for j in range(i+1):
                firstGroup[j] = -1
            for j in range(i, -1, -1):
                nextGroup[j] = firstGroup[inGroup[j]]
                firstGroup[inGroup[j]] = j
        #
        # Renumber to get rid of the numbers which were skipped
        #
        renumbered = np.zeros(nTargets, dtype='bool')
        nTmp = nGroups
        nGroups = 0
        for i in range(nTargets):
            if not renumbered[i]:
                j = firstGroup[inGroup[i]]
                while j != -1:
                    inGroup[j] = nGroups
                    renumbered[j] = True
                    j = nextGroup[j]
                nGroups += 1
        #
        # Reset the values of firstGroup and inGroup
        #
        firstGroup[:] = -1
        for i in range(nTargets-1, -1, -1):
            nextGroup[i] = firstGroup[inGroup[i]]
            firstGroup[inGroup[i]] = i
        #
        # Get the multiplicity
        #
        for i in range(nGroups):
            multGroup[i] = 0
            j = firstGroup[i]
            while j != -1:
                multGroup[i] += 1
                j = nextGroup[j]
        #
        # Set attributes
        #
        self.nGroups = nGroups
        self.nTargets = nTargets
        self.inGroup = inGroup
        self.multGroup = multGroup
        self.firstGroup = firstGroup
        self.nextGroup = nextGroup
        return


def spheregroup(ra, dec, linklength, chunksize=None):
    """Perform friends-of-friends grouping given ra/dec coordinates.

    Parameters
    ----------
    ra, dec : :class:`numpy.ndarray`
        Arrays of coordinates to group in decimal degrees.
    linklength : :class:`float`
        Linking length for the groups in decimal degrees.
    chunksize : :class:`float`, optional
        Break up the sphere into chunks of this size in decimal degrees.

    Returns
    -------
    :func:`tuple`
        A tuple containing the group number of each object, the multiplicity
        of each group, the first member of each group, and the next
        member of the group for each object.

    Raises
    ------
    :exc:`PydlutilsException`
        If the array of coordinates only contains one point.

    Notes
    -----
    It is important that `chunksize` >= 4 * `linklength`.  This is enforced.

    .. warning:: Behavior at the poles is not well tested.
    """
    npoints = ra.size
    if npoints == 1:
        raise PydlutilsException("Cannot group only one point!")
    #
    # Define the chunksize
    #
    if chunksize is not None:
        if chunksize < 4.0*linklength:
            chunksize = 4.0*linklength
            warn("chunksize changed to {0:.2f}.".format(chunksize), PydlutilsUserWarning)
    else:
        chunksize = max(4.0*linklength, 0.1)
    #
    # Initialize chunks
    #
    chunk = chunks(ra, dec, chunksize)
    chunk.assign(ra, dec, linklength)
    #
    # Run friends-of-friends
    #
    ingroup, multgroup, firstgroup, nextgroup, ngroups = chunk.friendsoffriends(ra, dec, linklength)
    #
    # Renumber the groups in order of appearance
    #
    renumbered = np.zeros(npoints, dtype='bool')
    iclump = 0
    for i in range(npoints):
        if not renumbered[i]:
            j = firstgroup[ingroup[i]]
            while j != -1:
                ingroup[j] = iclump
                renumbered[j] = True
                j = nextgroup[j]
            iclump += 1
    #
    # Reset the index lists
    #
    firstgroup[:] = -1
    for i in range(npoints-1, -1, -1):
        nextgroup[i] = firstgroup[ingroup[i]]
        firstgroup[ingroup[i]] = i
    #
    # Reset the multiplicities
    #
    multgroup[:] = 0
    for i in range(ngroups):
        j = firstgroup[i]
        while j != -1:
            multgroup[i] += 1
            j = nextgroup[j]
    return (ingroup, multgroup, firstgroup, nextgroup)


def spherematch(ra1, dec1, ra2, dec2, matchlength, chunksize=None,
                maxmatch=1):
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
    :func:`tuple`
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
    #
    # Set default values
    #
    if chunksize is None:
        chunksize = max(4.0*matchlength, 0.1)
    #
    # Check input size
    #
    if ra1.size == 1:
        raise PydlutilsException("Change the order of the sets of coordinates!")
    #
    # Initialize chunks
    #
    chunk = chunks(ra1, dec1, chunksize)
    chunk.assign(ra2, dec2, matchlength)
    #
    # Create return arrays
    #
    match1 = list()
    match2 = list()
    distance12 = list()
    for i in range(ra1.size):
        currra = np.fmod(ra1[i]+chunk.raOffset, 360.0)
        rachunk, decchunk = chunk.get(currra, dec1[i])
        jmax = len(chunk.chunkList[decchunk][rachunk])
        if jmax > 0:
            for j in range(jmax):
                k = chunk.chunkList[decchunk][rachunk][j]
                sep = gcirc(ra1[i], dec1[i], ra2[k], dec2[k], units=2)/3600.0
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
        gotten1 = np.zeros(ra1.size, dtype='i4')
        gotten2 = np.zeros(ra2.size, dtype='i4')
        nmatch = 0
        for i in range(omatch1.size):
            if (gotten1[omatch1[s[i]]] < maxmatch and
                    gotten2[omatch2[s[i]]] < maxmatch):
                gotten1[omatch1[s[i]]] += 1
                gotten2[omatch2[s[i]]] += 1
                nmatch += 1
        match1 = np.zeros(nmatch, dtype='i4')
        match2 = np.zeros(nmatch, dtype='i4')
        distance12 = np.zeros(nmatch, dtype='d')
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
    return (match1, match2, distance12)
