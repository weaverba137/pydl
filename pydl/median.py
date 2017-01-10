# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-


def median(array, width=None, axis=None, even=False):
    """Replicate the IDL ``MEDIAN()`` function.

    Parameters
    ----------
    array : array-like
        Compute the median of this array.
    width : :class:`int`, optional
        Size of the neighborhood in which to compute the median (*i.e.*,
        perform median filtering).  If omitted, the median of the whole
        array is returned.
    axis : :class:`int`, optional
        Compute the median over this axis for a multi-dimensional array.  If
        ommitted, the median over the entire array will be returned.  If
        set, this function will behave as though `even` is ``True``.
    even : :class:`bool`, optional
        If set to ``True``, the median of arrays with an even number of elements
        will be the average of the middle two values.

    Returns
    -------
    array-like
        The median of the array.

    Raises
    ------
    ValueError
        If `width` is set, and the input `array` is not 1 or 2 dimensional.

    Notes
    -----
    * For arrays with an even number of elements, the :func:`numpy.median`
      function behaves like ``MEDIAN(array, /EVEN)``, so the absence of
      the `even` keyword has to turn *off* that behavior.
    * For median filtering, this uses :func:`scipy.signal.medfilt` and
      :func:`scipy.signal.medfilt2d` under the hood, but patches up the
      values on the array boundaries to match the return values of the
      IDL ``MEDIAN()`` function.
    """
    import numpy as np
    from scipy.signal import medfilt, medfilt2d
    if width is None:
        if axis is None:
            f = array.flatten()
            if f.size % 2 == 1 or even:
                return np.median(array)
            else:
                i = f.argsort()
                return f[i[f.size//2]]
        else:
            return np.median(array, axis=axis)
    else:
        if array.ndim == 1:
            medarray = medfilt(array, min(width, array.size))
            istart = int((width - 1)/2)
            iend = array.size - int((width + 1)/2)
            i = np.arange(array.size)
            w = (i < istart) | (i > iend)
            medarray[w] = array[w]
            return medarray
        elif array.ndim == 2:
            medarray = medfilt2d(array, min(width, array.size))
            istart = int((width-1)/2)
            iend = (array.shape[0] - int((width+1)/2), array.shape[1] - int((width+1)/2))
            i = np.arange(array.shape[0])
            j = np.arange(array.shape[1])
            w = ((i < istart) | (i > iend[0]), (j < istart) | (j > iend[1]))
            medarray[w[0], :] = array[w[0], :]
            medarray[:, w[1]] = array[:, w[1]]
            return medarray
        else:
            raise ValueError("Invalid number of dimensions for input array!")
