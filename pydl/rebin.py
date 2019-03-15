# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-


def rebin(x, d, sample=False):
    """Resize `x` to new dimensions given by `d`.  The new dimensions must
    be integer multiples or factors of the original dimensions.

    Although there are some elegant solutions out there for rebinning, this
    function is intended to replace the IDL ``REBIN()`` function, which
    has a number of special properties:

    * It refuses to perform extrapolation when rebinning to a larger size
      in a particular dimension.
    * It can simultaneously rebin to a larger size in one dimension while
      rebinning to a smaller size in another dimension.

    Parameters
    ----------
    x : :class:`~numpy.ndarray`
        The array to resample.
    d : :func:`tuple`
        The new shape of the array.
    sample : :class:`bool`, optional
        If ``True``, nearest-neighbor techniques will be used instead of
        interpolation.

    Returns
    -------
    :class:`~numpy.ndarray`
        The resampled array.

    Raises
    ------
    :exc:`ValueError`
        If the new dimensions are incompatible with the algorithm.

    References
    ----------
    http://www.harrisgeospatial.com/docs/rebin.html

    Examples
    --------
    >>> from numpy import arange, float
    >>> from pydl import rebin
    >>> rebin(arange(10, dtype=float), (5,)) # doctest: +NORMALIZE_WHITESPACE
    array([ 0.5,  2.5,  4.5,  6.5,  8.5])
    >>> rebin(arange(5, dtype=float), (10,)) # doctest: +NORMALIZE_WHITESPACE
    array([ 0. ,  0.5,  1. ,  1.5,  2. ,  2.5,  3. ,  3.5,  4. ,  4. ])
    """
    from numpy import floor, zeros
    d0 = x.shape
    if len(d0) != len(d):
        raise ValueError(("The new shape is incompatible with the " +
                          "original array.!"))
    for k in range(len(d0)):
        if d[k] > d0[k]:
            if d[k] % d0[k] != 0:
                raise ValueError(("{0:d} is not a multiple " +
                                  "of {1:d}!").format(d[k], d0[k]))
        elif d[k] == d0[k]:
            pass
        else:
            if d0[k] % d[k] != 0:
                raise ValueError(("{0:d} is not a multiple " +
                                  "of {1:d}!").format(d0[k], d[k]))
    xx = x.copy()
    new_shape = list(d0)
    for k in range(len(d0)):
        new_shape[k] = d[k]
        r = zeros(new_shape, dtype=xx.dtype)
        sliceobj0 = [slice(None)]*len(d0)
        sliceobj1 = [slice(None)]*len(d0)
        sliceobj = [slice(None)]*len(d)
        f = d0[k]/d[k]
        if d[k] > d0[k]:
            for i in range(d[k]):
                p = f*i
                fp = int(floor(p))
                sliceobj0[k] = slice(fp, fp + 1)
                sliceobj[k] = slice(i, i + 1)
                if sample:
                    r[sliceobj] = xx[sliceobj0]
                else:
                    if p < d0[k] - 1:
                        sliceobj1[k] = slice(fp + 1, fp + 2)
                        rshape = r[sliceobj].shape
                        r[sliceobj] = (xx[sliceobj0].reshape(rshape) +
                                       (p - fp)*(xx[sliceobj1] -
                                                 xx[sliceobj0]).reshape(rshape)
                                      )
                    else:
                        r[sliceobj] = xx[sliceobj0]
        elif d[k] == d0[k]:
            for i in range(d[k]):
                sliceobj0[k] = slice(i, i + 1)
                sliceobj[k] = slice(i, i + 1)
                r[sliceobj] = xx[sliceobj0]
        else:
            for i in range(d[k]):
                sliceobj[k] = slice(i, i + 1)
                if sample:
                    fp = int(floor(f*i))
                    sliceobj0[k] = slice(fp, fp + 1)
                    r[sliceobj] = xx[sliceobj0]
                else:
                    sliceobj0[k] = slice(int(f*i), int(f*(i+1)))
                    rshape = r[sliceobj].shape
                    r[sliceobj] = xx[sliceobj0].sum(k).reshape(rshape)/f
        xx = r
    return r
