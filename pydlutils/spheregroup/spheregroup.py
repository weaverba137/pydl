#
# $Id$
#
def spheregroup(ra,dec,linklength,**kwargs):
    """Perform friends-of-friends grouping given ra/dec coordinates
    """
    npoints = ra.size
    #
    # Define the chunksize
    #
    if 'chunksize' in kwargs:
        if kwargs['chunksize'] < 4.0*linklength:
            chunksize = 4.0*linklength
            print "chunksize changed to %f." % chunksize
        else:
            chunksize = kwargs['chunksize']
    else:
        chunksize = max(4.0*linklength,0.1)
    #
    # Initialize chunks
    #
    chunk = chunks(ra,dec,chunksize)
    chunk.assign(ra,dec,linklength)
    if 'debug' in kwargs:
        print "raOffset = %7.3f." % chunk.raOffset
        print chunk.raBounds
        print chunk.decBounds
        print chunk.chunkList
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
