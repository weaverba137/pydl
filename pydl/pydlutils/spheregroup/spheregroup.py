# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from __future__ import print_function
#
def spheregroup(ra,dec,linklength,chunksize=None):
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
    spheregroup : :func:`tuple`
        A tuple containing the group number of each object, the multiplicity
        of each group, the first member of each group, and the next
        member of the group for each object.

    Raises
    ------
    PydlutilsException
        If the array of coordinates only contains one point.

    Notes
    -----
    It is important that `chunksize` >= 4 * `linklength`.  This is enforced.

    .. warning:: Behavior at the poles is not well tested.
    """
    from . import chunks
    from .. import PydlutilsException
    from .. import PydlutilsUserWarning
    from numpy import zeros
    from warnings import warn
    npoints = ra.size
    if npoints == 1:
        raise PydlutilsException("Cannot group only one point!")
    #
    # Define the chunksize
    #
    if chunksize is not None:
        if chunksize < 4.0*linklength:
            chunksize = 4.0*linklength
            warn("chunksize changed to {0:.2f}.".format(chunksize),PydlutilsUserWarning)
    else:
        chunksize = max(4.0*linklength,0.1)
    #
    # Initialize chunks
    #
    chunk = chunks(ra,dec,chunksize)
    chunk.assign(ra,dec,linklength)
    #
    # Run friends-of-friends
    #
    ingroup,multgroup,firstgroup,nextgroup,ngroups = chunk.friendsoffriends(ra,dec,linklength)
    #
    # Renumber the groups in order of appearance
    #
    renumbered = zeros(npoints,dtype='bool')
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
    multgroup[:] = 0
    for i in range(ngroups):
        j = firstgroup[i]
        while j != -1:
            multgroup[i] += 1
            j = nextgroup[j]
    return (ingroup,multgroup,firstgroup,nextgroup)
