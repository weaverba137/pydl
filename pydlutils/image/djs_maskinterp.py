#
#
#
def djs_maskinterp1(yval,mask,xval=None,**kwargs):
    """Helper program for djs_maskinterp.
    """
    import numpy as np
    good = mask == 0
    if good.all():
        return yval
    ngood = good.sum()
    igood = good.nonzero()[0]
    if ngood == 0:
        return yval
    if ngood == 1:
        return np.zeros(yval.shape,dtype=yval.dtype) + yval[igood[0]]
    ynew = yval.astype('d')
    ny = yval.size
    ibad = (mask != 0).nonzero()[0]
    if xval is None:
        ynew[ibad] = np.interp(ibad,igood,ynew[igood])
        if 'const' in kwargs:
            if igood[0] != 0:
                ynew[0:igood[0]] = ynew[igood[0]]
            if igood[ngood-1] != ny-1:
                ynew[igood[ngood-1]+1:ny] = ynew[igood[ngood-1]]
    else:
        ii = xval.argsort()
        ibad = (mask[ii] != 0).nonzero()[0]
        igood = (mask[ii] == 0).nonzero()[0]
        ynew[ii[ibad]] = np.interp(xval[ii[ibad]],xval[ii[igood]],ynew[ii[igood]])
        if 'const' in kwargs:
            if igood[0] != 0:
                ynew[ii[0:igood[0]]] = ynew[ii[igood[0]]]
            if igood[ngood-1] != ny-1:
                ynew[ii[igood[ngood-1]+1:ny]] = ynew[ii[igood[ngood-1]]]
    return ynew

#
#
#
def djs_maskinterp(yval,mask,xval=None,**kwargs):
    """Interpolate over masked pixels in a vector, image or 3-D array.

    Arguments:
    yval -- The input values
    mask -- The mask
    xval -- If set, use these x values, otherwise use an array
    """
    import numpy as np
    if mask.shape != yval.shape:
        raise ValueError('mask must have the same shape as yval.')
    if xval is not None:
        if xval.shape != yval.shape:
            raise ValueError('xval must have the same shape as yval.')
    ndim = yval.ndim
    if ndim == 1:
        ynew = djs_maskinterp1(yval,mask,xval,**kwargs)
    else:
        if 'axis' not in kwargs:
            raise ValueError('Must set axis if yval has more than one dimension.')
        axis = kwargs['axis']
        if axis < 0 or axis > ndim-1 or axis - int(axis) != 0:
            raise ValueError('Invalid axis value.')
        ynew = np.zeros(yval.shape,dtype=yval.dtype)
        if ndim == 2:
            if xval is None:
                if axis == 0:
                    for i in range(yval.shape[0]):
                        ynew[i,:] = djs_maskinterp1(yval[i,:],mask[i,:],**kwargs)
                else:
                    for i in range(yval.shape[1]):
                        ynew[:,i] = djs_maskinterp1(yval[:,i],mask[:,i],**kwargs)
            else:
                if axis == 0:
                    for i in range(yval.shape[0]):
                        ynew[i,:] = djs_maskinterp1(yval[i,:],mask[i,:],xval[i,:],**kwargs)
                else:
                    for i in range(yval.shape[1]):
                        ynew[:,i] = djs_maskinterp1(yval[:,i],mask[:,i],xval[:,i],**kwargs)
        elif ndim == 3:
            if xval is None:
                if axis == 0:
                    for i in range(yval.shape[0]):
                        for j in range(yval.shape[1]):
                            ynew[i,j,:] = djs_maskinterp1(yval[i,j,:],mask[i,j,:],**kwargs)
                elif axis == 1:
                    for i in range(yval.shape[0]):
                        for j in range(yval.shape[2]):
                            ynew[i,:,j] = djs_maskinterp1(yval[i,:,j],mask[i,:,j],**kwargs)
                else:
                    for i in range(yval.shape[1]):
                        for j in range(yval.shape[2]):
                            ynew[:,i,j] = djs_maskinterp1(yval[:,i,j],mask[:,i,j],**kwargs)
            else:
                if axis == 0:
                    for i in range(yval.shape[0]):
                        for j in range(yval.shape[1]):
                            ynew[i,j,:] = djs_maskinterp1(yval[i,j,:],mask[i,j,:],xval[i,j,:],**kwargs)
                elif axis == 1:
                    for i in range(yval.shape[0]):
                        for j in range(yval.shape[2]):
                            ynew[i,:,j] = djs_maskinterp1(yval[i,:,j],mask[i,:,j],xval[i,j,:],**kwargs)
                else:
                    for i in range(yval.shape[1]):
                        for j in range(yval.shape[2]):
                            ynew[:,i,j] = djs_maskinterp1(yval[:,i,j],mask[:,i,j],xval[i,j,:],**kwargs)
        else:
            raise NotImplementedError('Unsupported number of dimensions.')
    return ynew


