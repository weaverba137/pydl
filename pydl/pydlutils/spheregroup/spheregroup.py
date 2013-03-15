# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from __future__ import print_function
#
def spheregroup(ra,dec,linklength,chunksize=None,debug=False):
    """Perform friends-of-friends grouping given ra/dec coordinates

    Parameters
    ----------
    ra, dec : ndarray
        Arrays of coordinates to group.
    linklength : float
    chunksize : float
    debug : bool

    Returns
    -------
    spheregroup : tuple
    """
    npoints = ra.size
    #
    # Define the chunksize
    #
    if chunksize is not None:
        if chunksize < 4.0*linklength:
            chunksize = 4.0*linklength
            print("chunksize changed to {0:.2f}.".format(chunksize))
    else:
        chunksize = max(4.0*linklength,0.1)
    #
    # Initialize chunks
    #
    chunk = chunks(ra,dec,chunksize)
    chunk.assign(ra,dec,linklength)
    if debug:
        print("raOffset = {0:7.3f}.".format(chunk.raOffset))
        print(chunk.raBounds)
        print(chunk.decBounds)
        print(chunk.chunkList)
    #
    # Run friends-of-friends
    #
    ingroup,multgroup,firstgroup,nextgroup,ngroups = chunk.friendsoffriends(ra,dec,linklength)
    #
    # Renumber the groups in order of appearance
    #
    renumbered = self.np.zeros(npoints,dtype='bool')
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
    for i in range(npoints-1,-1,-1):
        nextgroup[i] = firstgroup[ingroup[i]]
        firstgroup[ingroup[i]] = i
    #
    # Reset the multiplicities
    #
    for i in range(ngroups):
        j = firstgroup[i]
        while j != -1:
            multgroup[i] += 1
            j = nextgroup[j]
    return (ingroup,multgroup,firstgroup,nextgroup)
