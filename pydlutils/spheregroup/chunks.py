#
# $Id$
#
class chunks(object):
    """chunks class

    Functions for creating and manipulating spherical chunks are implemented
    as methods on this class.
    """
    import numpy as np
    from pydlutils import PydlutilsException
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
        self.nDec = 3 + int(self.np.floor(decRange/minSize))
        decRange = minSize*float(self.nDec)
        decMin = decMin - 0.5*(decRange - decMax + decMin)
        decMax = decMin + decRange
        if decMin < -90.0 + 3.0*minSize:
            decMin = -90.0
        if decMax > 90.0 - 3.0*minSize:
            decMax = 90.0
        self.decBounds = decMin + (decMax-decMin)*self.np.arange(self.nDec+1,dtype='d')/float(self.nDec)
        #
        # Find ra offset which minimizes the range in ra (this should take care
        # of the case that ra crosses zero in some parts
        #
        if abs(self.decBounds[self.nDec]) > abs(self.decBounds[0]):
            cosDecMin = self.np.cos(self.np.deg2rad(self.decBounds[self.nDec]))
        else:
            cosDecMin = self.np.cos(self.np.deg2rad(self.decBounds[0]))
        if cosDecMin <= 0.0:
            raise self.PydlutilsException("cosDecMin=%f not positive in setchunks()." % cosDecMin)
        self.raRange, self.raOffset = self.rarange(ra,minSize/cosDecMin)
        self.raMin, self.raMax = self.getraminmax(ra,self.raOffset)
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
                cosDecMin = self.np.cos(self.np.deg2rad(self.decBounds[i]))
            else:
                cosDecMin = self.np.cos(self.np.deg2rad(self.decBounds[i+1]))
            if cosDecMin <= 0.0:
                raise PydlutilsException("cosDecMin=%f not positive in setchunks()." % cosDecMin)
            #
            # Get raBounds array for this declination array, leave an extra
            # cell on each end
            #
            self.nRa.append(3 + int(self.np.floor(cosDecMin*self.raRange/minSize)))
            raRangeTmp = minSize*float(self.nRa[i])/cosDecMin
            raMinTmp = self.raMin - 0.5*(raRangeTmp-self.raMax+self.raMin)
            raMaxTmp = raMinTmp + raRangeTmp
            #
            # If we cannot avoid the 0/360 point, embrace it
            #
            if ( raRangeTmp >= 360.0 or
                raMinTmp <= minSize/cosDecMin or
                raMaxTmp >= 360.0 - minSize/cosDecMin or
                abs(self.decBounds[i]) == 90.0 ):
                raMinTmp = 0.0
                raMaxTmp = 360.0
                raRangeTmp = 360.0
            if self.decBounds[i] == -90.0 or self.decBounds[i+1] == 90.0:
                self.nRa[i] = 1
            self.raBounds.append(raMinTmp +
                (raMaxTmp-raMinTmp)*self.np.arange(self.nRa[i]+1,dtype='d')/
                float(self.nRa[i]))
        #
        # Create an empty set of lists to hold the output of self.assign()
        #
        self.chunkList = [[list() for j in range(self.nRa[i])] for i in range(self.nDec)]
        return
    def rarange(self,ra,minSize):
        """Finds the offset which yields the smallest raRange & returns both.
        """
        NRA = 6
        raRangeMin = 361.0
        raOffset = 0.0
        EPS = 1.0e-5
        for j in range(NRA):
            raMin,raMax = self.getraminmax(ra, 360.0*float(j)/float(NRA))
            raRange = raMax-raMin
            if (2.0*(raRange-raRangeMin)/(raRange+raRangeMin) < -EPS and
                raMin > minSize and raMax < 360.0 - minSize):
                raRangeMin = raRange
                raOffset = 360.0*float(j)/float(NRA)
        return (raRangeMin, raOffset)
    def getraminmax(self,ra,raOffset):
        """Utility function used by rarange.
        """
        currRa = self.np.fmod(ra + raOffset, 360.0)
        return (currRa.min(), currRa.max())
    def cosDecMin(self,i):
        """Frequently used utility function.
        """
        if abs(self.decBounds[i]) > abs(self.decBounds[i+1]):
            return self.np.cos(self.np.deg2rad(self.decBounds[i]))
        else:
            return self.np.cos(self.np.deg2rad(self.decBounds[i+1]))
    def assign(self, ra, dec, marginSize):
        """Take the objects and the chunks (already defined in the constructor)
        and assign the objects to the appropriate chunks, with some leeway
        given by the parameter marginSize.  Basically, at the end, each
        chunk should be associated with a list of the objects that belong
        to it.
        """
        if marginSize >= self.minSize:
            raise self.PydlutilsException("marginSize>=minSize (%f>=%f) in chunks.assign()." % (marginSize,self.minSize))
        chunkDone = [[False for j in range(self.nRa[i])] for i in range(self.nDec)]
        for i in range(ra.size):
            currRa = self.np.fmod(ra[i]+self.raOffset,360.0)
            try:
                raChunkMin, raChunkMax, decChunkMin, decChunkMax = self.getbounds(currRa,dec[i],marginSize)
            except self.PydlutilsException:
                continue
            #
            # Reset chunkDone.  This is silly, but is necessary to
            # reproduce the logic.
            #
            for decChunk in range(decChunkMin,decChunkMax+1):
                for raChunk in range(raChunkMin[decChunk-decChunkMin]-1,raChunkMax[decChunk-decChunkMin]+2):
                    if raChunk < 0:
                        currRaChunk = (raChunk+self.nRa[decChunk]) % self.nRa[decChunk]
                    elif raChunk > self.nRa[decChunk]-1:
                        currRaChunk = (raChunk-self.nRa[decChunk]) % self.nRa[decChunk]
                    else:
                        currRaChunk = raChunk
                    if currRaChunk >= 0 and currRaChunk <= self.nRa[decChunk]-1:
                        chunkDone[decChunk][currRaChunk] = False
            for decChunk in range(decChunkMin,decChunkMax+1):
                for raChunk in range(raChunkMin[decChunk-decChunkMin],raChunkMax[decChunk-decChunkMin]+1):
                    if raChunk < 0:
                        currRaChunk = (raChunk+self.nRa[decChunk]) % self.nRa[decChunk]
                    elif raChunk > self.nRa[decChunk]-1:
                        currRaChunk = (raChunk-self.nRa[decChunk]) % self.nRa[decChunk]
                    else:
                        currRaChunk = raChunk
                    if currRaChunk >= 0 and currRaChunk <= self.nRa[decChunk]-1:
                        if not chunkDone[decChunk][currRaChunk]:
                            self.chunkList[decChunk][currRaChunk].append(i)
                            chunkDone[decChunk][currRaChunk] = True
        return
    def getbounds(self,ra,dec,marginSize):
        """Find the set of chunks a point (with margin) belongs to.
        """
        #
        # Find the declination slice without regard to marginSize
        #
        decChunkMin = int(self.np.floor((dec - self.decBounds[0])*
            float(self.nDec)/
            (self.decBounds[self.nDec]-self.decBounds[0])))
        decChunkMax = decChunkMin
        if decChunkMin < 0 or decChunkMin > self.nDec - 1:
            raise self.PydlutilsException("decChunkMin out of range in chunks.getbounds().")
        #
        # Set minimum and maximum bounds of dec
        #
        while dec-self.decBounds[decChunkMin] < marginSize and decChunkMin > 0:
            decChunkMin -= 1
        while self.decBounds[decChunkMax+1]-dec < marginSize and decChunkMax < self.nDec -1:
            decChunkMax += 1
        #
        # Find ra chunk bounds for each dec chunk
        #
        raChunkMin = self.np.zeros(decChunkMax-decChunkMin+1,dtype='i4')
        raChunkMax = self.np.zeros(decChunkMax-decChunkMin+1,dtype='i4')
        for i in range(decChunkMin,decChunkMax+1):
            cosDecMin = self.cosDecMin(i)
            raChunkMin[i-decChunkMin] = int(self.np.floor((ra - self.raBounds[i][0])*
                float(self.nRa[i])/
                (self.raBounds[i][self.nRa[i]] - self.raBounds[i][0])))
            raChunkMax[i-decChunkMin] = raChunkMin[i-decChunkMin]
            if raChunkMin[i-decChunkMin] < 0 or raChunkMin[i-decChunkMin] > self.nRa[i]-1:
                raise self.PydlutilsException("raChunkMin out of range in chunks.getbounds().")
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
        decChunk = int(self.np.floor((dec - self.decBounds[0])*
            float(self.nDec)/
            (self.decBounds[self.nDec]-self.decBounds[0])))
        #
        # Find ra chunk
        #
        if decChunk < self.nDec and decChunk >= 0:
            raChunk = int(self.np.floor((ra - self.raBounds[decChunk][0])*
                float(self.nRa[decChunk])/
                (self.raBounds[decChunk][self.nRa[decChunk]] - self.raBounds[decChunk][0])))
            if raChunk < 0 or raChunk > self.nRa[decChunk]-1:
                raise self.PydlutilsException("raChunk out of range in chunks.get()")
        else:
            raChunk = -1
        return (raChunk, decChunk)
