# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def djs_median(array,dimension=None,width=None,boundary='none'):
    """Compute the median of an array.

    Use a filtering box or collapse the image along one dimension.

    Parameters
    ----------
    array : ndarray
        input array
    dimension : int, optional
        Compute the median over this dimension. It is an error to specify both
        `dimension` and `width`.
    width : int, optional
        Width of the median window. It is an error to specify both
        `dimension` and `width`.
    boundary : { 'none', 'reflect', 'nearest', 'wrap' }, optional
        Boundary condition to impose.  'none' means no filtering is done within
        `width`/2 of the boundary.  'reflect' means reflect pixel values around the
        boundary. 'nearest' means use the values of the nearest boundary pixel.
        'wrap' means wrap pixel values around the boundary. 'nearest' and 'wrap'
        are not implemented.

    Returns
    -------
    djs_median : ndarray
        The output.  If neither `dimension` nor `width` are set, this is a scalar
        value, just the output of ``numpy.median()``.  If `width` is set, the
        result has the same shape as the input array.
    """
    import numpy as np
    from scipy.signal import medfilt
    if dimension is None and width is None:
        return np.median(array)
    elif dimension is None:
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
    elif width is None:
        raise NotImplementedError('This type of median is not implemented')
    else:
        raise ValueError('Invalid to specify both dimension & width.')
    return medarray
