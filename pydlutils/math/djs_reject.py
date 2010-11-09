#
#
#
def djs_reject(data,model,**kwargs):
    """Routine to reject points when doing an iterative fit to data.

    Arguments:
    data -- The data
    model -- The model
    """
    import numpy as np
    from pydlutils.misc import djs_laxisnum
    if model is not None:
        if data.shape != model.shape:
            raise ValueError('Dimensions of data and model do not agree.')
    if 'inmask' in kwargs:
        if data.shape != kwargs['inmask'].shape:
            raise ValueError('Dimensions of data and inmask do not agree.')
    if 'maxrej' in kwargs:
        if 'groupdim' in kwargs:
            if len(kwargs['maxreg']) != len(kwargs['groupdim']):
                raise ValueError('maxrej and groupdim must have the same number of elements.')
        if 'groupsize' in kwargs:
            if len(kwargs['maxreg']) != len(kwargs['groupsize']):
                raise ValueError('maxrej and groupsize must have the same number of elements.')
            groupsize = kwargs['groupsize']
        else:
            groupsize = len(data)
    #
    # Create outmask setting = 1 for good data.
    #
    if 'outmask' in kwargs:
        if kwargs['outmask'] is None:
            outmask = np.ones(data.shape,dtype='bool')
        else:
            if kwargs['outmask'].shape != data.shape:
                raise ValueError('Dimensions of data and outmask do not agree.')
            outmask = kwargs['outmask']
    else:
        raise ValueError('Try passing outmask=None.')
    if model is None:
        if 'inmask' in kwargs:
            outmask = kwargs['inmask']
        return (outmask,0)
    if 'sigma' in kwargs and 'invvar' in kwargs:
        raise ValueError('Cannot set both sigma and invvar.')
    if 'sigma' not in kwargs and 'invvar' not in kwargs:
        if 'inmask' in kwargs:
            igood = (kwargs['inmask'] & outmask).nonzero()[0]
        else:
            igood = outmask.nonzero()[0]
        if len(igood > 1):
            sigma = np.std(data[igood] - model[igood])
        else:
            sigma = 0
    diff = data - model
    #
    # The working array is badness, which is set to zero for good points
    # (or points already rejected), and positive values for bad points.
    # The values determine just how bad a point is, either corresponding
    # to the number of sigma above or below the fit, or to the number
    # of multiples of maxdev away from the fit.
    #
    badness = np.zeros(outmask.shape,dtype=data.dtype)
    #
    # Decide how bad a point is according to lower.
    #
    if 'lower' in kwargs:
        if 'invvar' in kwargs:
            qbad = (diff * np.sqrt(kwargs['invvar'])) < -kwargs['lower']
            badness += ((-diff * np.sqrt(kwargs['invvar'])) > 0) * qbad
        else:
            qbad = diff < (-lower * sigma)
            badness += ((-diff/(sigma + (sigma == 0))) > 0) * qbad
    #
    # Decide how bad a point is according to upper.
    #
    if 'upper' in kwargs:
        if 'invvar' in kwargs:
            qbad = (diff * np.sqrt(kwargs['invvar'])) > kwargs['upper']
            badness += ((diff * np.sqrt(kwargs['invvar'])) > 0) * qbad
        else:
            qbad = diff > (lower * sigma)
            badness += ((diff/(sigma + (sigma == 0))) > 0) * qbad
    #
    # Decide how bad a point is according to maxdev.
    #
    if 'maxdev' in kwargs:
        qbad = np.absolute(diff) > kwargs['maxdev']
        badness += np.absolute(diff) / kwargs['maxdev'] * qbad
    #
    # Do not consider rejecting points that are already rejected by inmask.
    # Do not consider rejecting points that are already rejected by outmask,
    # if sticky is set.
    #
    if 'inmask' in kwargs:
        badness *= kwargs['inmask']
    if 'sticky' in kwargs:
        badness *= outmask
    #
    # Reject a maximum of maxrej (additional) points in all the data, or
    # in each group as specified by groupsize, and optionally along each
    # dimension specified by groupdim.
    #
    if 'maxrej' in kwargs:
        #
        # Loop over each dimension of groupdim or loop once if not set.
        #
        for iloop in range(max(len(kwargs['groupdim']),1)):
            #
            # Assign an index number in this dimension to each data point.
            #
            if len(kwargs['groupdim']) > 0:
                yndim = len(ydata.shape)
                if kwargs['groupdim'][iloop] > yndim:
                    raise ValueError('groupdim is larger than the number of dimensions for ydata.')
                dimnum = djs_laxisnum(ydata.shape,iaxis=kwargs['groupdim'][iloop]-1)
            else:
                dimnum = 0
            #
            # Loop over each vector specified by groupdim. For example, if
            # this is a 2-D array with groupdim=1, then loop over each
            # column of the data.  If groupdim=2, then loop over each row.
            # If groupdim is not set, then use the whole image.
            #
            for ivec in range(max(dimnum)):
                #
                # At this point it is not possible that dimnum is not set.
                #
                indx = (dimnum == ivec).nonzero()[0]
                #
                # Within this group of points, break it down into groups
                # of points specified by groupsize, if set.
                #
                nin = len(indx)
                if 'groupbadpix' in kwargs:
                    goodtemp = badness == 0
                    groups_lower = (-1*np.diff(np.insert(goodtemp,0,1)) == 1).nonzero()[0]
                    groups_upper = (np.diff(np.append(goodtemp,1)) == 1).nonzero()[0]
                    ngroups = len(groups_lower)
                else:
                    if 'groupsize' in kwargs:
                        ngroups = nin/groupsize + 1
                        groups_lower = np.arange(ngroups,dtype='i4')*groupsize
                        foo = (np.arange(ngroups,dtype='i4')+1)*groupsize
                        groups_upper = np.where(foo < nin,foo,nin) -1
                    else:
                        ngroups = 1
                        groups_lower = [0,]
                        groups_upper = [nin - 1,]
                for igroup in range(ngroups):
                    i1 = groups_lower[igroup]
                    i2 = groups_upper[igroup]
                    nii = i2 - i1 + 1
                    #
                    # Need the test that i1 != -1 below to prevent a crash
                    # condition, but why is it that we ever get groups
                    # without any points?  Because this is badly-written,
                    # that's why.
                    #
                    if nii > 0 and i1 != -1:
                        jj = indx[i1:i2+1]
                        #
                        # Test if too many points rejected in this group.
                        #
                        if np.sum(badness[jj] != 0) > kwargs['maxrej'][iloop]:
                            isort = badness[jj].argsort()
                            #
                            # Make the following points good again.
                            #
                            badness[jj[isort[0:nii-kwargs['maxrej'][iloop]]]] = 0
                        i1 += groupsize[iloop]
    #
    # Now modify outmask, rejecting points specified by inmask=0, outmask=0
    # if sticky is set, or badness > 0.
    #
    # print badness
    newmask = badness == 0
    # print newmask
    if 'grow' in kwargs:
        rejects = newmask==0
        if rejects.any():
            irejects = rejects.nonzero()[0]
            for k in range(1,kwargs['grow']):
                newmask[(irejects - k) > 0] = 0
                newmask[(irejects + k) < (data.shape[0]-1)] = 0
    if 'inmask' in kwargs:
        newmask = newmask & kwargs['inmask']
    if 'sticky' in kwargs:
        newmask = newmask & outmask
    #
    # Set qdone if the input outmask is identical to the output outmask.
    #
    qdone = np.all(newmask == outmask)
    outmask = newmask
    return (outmask,qdone)

