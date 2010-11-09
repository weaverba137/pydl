#
#
#
def djs_laxisnum(dims,iaxis=0):
    import numpy as np
    ndimen = len(dims)
    result = np.zeros(dims,dtype='i4')
    if ndimen == 1:
        pass
    elif ndimen == 2:
        if iaxis == 0:
            for k in range(dims[0]):
                result[k,:] = k
        elif iaxis == 1:
            for k in range(dims[1]):
                result[:,k] = k
        else:
            raise ValueError("Bad value for iaxis: %d" % iaxis)
    elif ndimen == 3:
        if iaxis == 0:
            for k in range(dims[0]):
                result[k,:,:] = k
        elif iaxis == 1:
            for k in range(dims[1]):
                result[:,k,:] = k
        elif iaxis == 2:
            for k in range(dims[2]):
                result[:,:,k] = k
        else:
            raise ValueError("Bad value for iaxis: %d" % iaxis)
    else:
        raise NotImplementedError("%d dimensions not supported." % ndimen)
    return result

