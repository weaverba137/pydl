# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def find_contiguous(x):
    """Find the longest sequence of contiguous non-zero array elements.

    Parameters
    ----------
    x : :class:`numpy.ndarray`
        A 1d array. A dtype of bool is preferred although any dtype where the
        operation ``if x[k]:`` is well-defined should work.

    Returns
    -------
    find_contiguous : :class:`list`
        A list of indices of the longest contiguous non-zero sequence.

    Examples
    --------
    >>> import numpy as np
    >>> from pydl.pydlutils.math import find_contiguous
    >>> find_contiguous(np.array([0,1,1,1,0,1,1,0,1]))
    [1, 2, 3]
    """
    contig = list()
    for k in range(x.size):
        if x[k]:
            if len(contig) == 0:
                contig.append([k])
            else:
                if k == contig[-1][-1]+1:
                    contig[-1].append(k)
                else:
                    contig.append([k])
    lengths = [len(c) for c in contig]
    longest = contig[lengths.index(max(lengths))]
    return longest
