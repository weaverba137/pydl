# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-


def uniq(x, index=None):
    """Replicates the IDL ``UNIQ()`` function.

    Returns the *subscripts* of the unique elements of an array.  The elements
    must actually be *sorted* before being passed to this function.  This can
    be done by sorting `x` explicitly or by passing the array subscripts that
    sort `x` as a second parameter.

    Parameters
    ----------
    x : array-like
        Search this array for unique items.
    index : array-like, optional
        This array provides the array subscripts that sort `x`.

    Returns
    -------
    array-like
        The subscripts of `x` that are the unique elements of `x`.

    Notes
    -----
    Given a sorted array, and assuming that there is a set of
    adjacent identical items, ``uniq()`` will return the subscript of the
    *last* unique item.  This charming feature is retained for
    reproducibility.

    References
    ----------
    http://www.exelisvis.com/docs/UNIQ.html

    Examples
    --------
    >>> import numpy as np
    >>> from pydl import uniq
    >>> data = np.array([ 1, 2, 3, 1, 5, 6, 1, 7, 3, 2, 5, 9, 11, 1 ])
    >>> print(uniq(np.sort(data)))
    [ 3  5  7  9 10 11 12 13]
    """
    from numpy import array, roll
    if index is None:
        indicies = (x != roll(x, -1)).nonzero()[0]
        if indicies.size > 0:
            return indicies
        else:
            return array([x.size - 1, ])
    else:
        q = x[index]
        indicies = (q != roll(q, -1)).nonzero()[0]
        if indicies.size > 0:
            return index[indicies]
        else:
            return array([q.size - 1, ], dtype=index.dtype)
