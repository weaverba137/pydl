# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
#
def decode_mixed(x):
    """Convert bytes in Numpy arrays into strings.  Leave other stuff alone.

    Parameters
    ----------
    x : object
        Input object.

    Returns
    -------
    decode_mixed : object
        If ``x`` has a ``decode()`` method, ``x.decode()`` will be returned.
        Otherwise ``x`` will be returned unchanged.
    """
    try:
        return x.decode()
    except:
        return x
