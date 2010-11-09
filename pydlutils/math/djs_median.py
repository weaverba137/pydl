#
#
#
def djs_median(array,**kwargs):
    """Compute the median of an array.

    Arguments:
    array -- input array
    Keyword Arguments:
    width -- Do a median filter with this window size.
    axis  -- Perform the median along this axis.
    """
    import numpy as np
    from scipy.signal import medfilt
    if 'boundary' in kwargs:
        boundary = kwargs['boundary']
    else:
        boundary = 'none'
    if 'axis' not in kwargs and 'width' not in kwargs:
        return np.median(array)
    elif 'axis' not in kwargs:
        width = kwargs['width']
        if boundary == 'none':
            if width == 1:
                medarray = array
            else:
                medarray = medfilt(array,min(width,array.size))
                #
                # medfilt doesn't handle edges the same way as IDL, so we
                # have to fix them up
                #
                istart = (width-1)/2
                iend = array.size - (width+1)/2
                i = np.arange(array.size)
                w = (i < istart) | (i > iend)
                medarray[w] = array[w]
        else:
            padsize = np.ceil(width/2.0)
            if array.ndim == 1:
                bigarr = np.zeros(array.shape[0]+2*padsize,dtype=array.dtype)
                bigarr[padsize:padsize+array.shape[0]] = array
            elif array.ndim == 2:
                bigarr = np.zeros((array.shape[0]+2*padsize,array.shape[1]+2*padsize),dtype=array.dtype)
                bigarr[padsize:padsize+array.shape[0],padsize:padsize+array.shape[1]] = array
            else:
                raise ValueError('Unsupported number of dimensions with this boundary condition.')
            if array.ndim == 1:
                if width == 1:
                    medarray = array
                else:
                    bigarr[0:padsize] = array[0:padsize][::-1]
                    bigarr[padsize+array.shape[0]:padsize*2+array.shape[0]] = array[array.shape[0]-padsize:array.shape[0]]
                    f = medfilt(bigarr,width)
                    istart = (width-1)/2
                    iend = array.size - (width+1)/2
                    i = np.arange(array.size)
                    w = (i < istart) | (i > iend)
                    f[w] = array[w]
                    medarray = f[padsize:padsize+array.shape[0]]
            else:
                if boundary == 'nearest':
                    raise NotImplementedError('This boundary condition not implemented')
                elif boundary == 'wrap':
                    raise NotImplementedError('This boundary condition not implemented')
                elif boundary == 'reflect':
                    # medarray = array
                    raise NotImplementedError('This boundary condition not implemented')
                else:
                    raise ValueError('Unknown boundary condition.')
    elif 'width' not in kwargs:
        raise NotImplementedError('This type of median is not implemented')
    else:
        raise ValueError('Invalid to specify both axis & width.')
    return medarray

