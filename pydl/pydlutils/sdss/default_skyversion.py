# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def default_skyversion():
    """Returns skyversion number to use for photoop imaging.

    Returns
    -------
    default_skyversion : :class:`int`
        The default skyversion number.

    Notes
    -----
    The skyversion number is hard-coded to 2.

    Examples
    --------
    >>> from pydl.pydlutils.sdss import default_skyversion
    >>> default_skyversion()
    2
    """
    return 2
