# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from __future__ import print_function
#
class groups(object):
    """Group a set of objects (a list of coordinates in some space) based on
    a friends-of-friends algorithm
    """
    import numpy as np
    from .. import PydlutilsException
    @staticmethod
    def euclid(x1,x2):
        """Pythagorean theorem in Euclidean space with arbitrary number
        of dimensions.
        """
        return self.np.sqrt(((x1-x2)**2).sum())
    @staticmethod
    def sphereradec(x1,x2):
        """Separation of two points on a 2D-sphere, assuming they are in
        longitude-latitude or right ascension-declination form.  Assumes
        everything is already in radians.
        """
        from ...goddard.astro import gcirc
        return gcirc(x1[0],x1[1],x2[0],x2[1],units=0)
    def __init__(self,coordinates,distance,separation='euclid'):
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
        elif isinstance(separation,str):
            if separation == 'euclid':
                self.separation = self.euclid
            elif separation == 'sphereradec':
                self.separation = self.sphereradec
            else:
                raise self.PydlutilsException("Unknown separation function: {0}.".format(separation))
        else:
            raise self.PydlutilsException("Improper type for separation!")
        #
        # Save information about the coordinates.
        #
        nGroups = 0
        nTargets = coordinates.shape[1]
        multGroup = self.np.zeros(nTargets,dtype='i4')
        firstGroup = self.np.zeros(nTargets,dtype='i4') -1
        nextGroup = self.np.zeros(nTargets,dtype='i4') -1
        inGroup = self.np.arange(nTargets,dtype='i4')
        #
        # Find all the other targets associated with each target
        #
        for i in range(nTargets):
            nTmp = 0
            minGroup = nGroups
            for j in range(nTargets):
                sep = self.separation(coordinates[:,i],coordinates[:,j])
                if sep <= distance:
                    multGroup[nTmp] = j
                    minGroup = min(minGroup,inGroup[j])
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
            for j in range(i,-1,-1):
                nextGroup[j] = firstGroup[inGroup[j]]
                firstGroup[inGroup[j]] = j
        #
        # Renumber to get rid of the numbers which were skipped
        #
        renumbered = self.np.zeros(nTargets,dtype='bool')
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
        for i in range(nTargets-1,-1,-1):
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
