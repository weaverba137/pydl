# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
====
pydl
====

Python replacements for functions that are part of the IDL_ built-in library, or
part of astronomical IDL libraries.  The emphasis is on reproducing results of
the astronomical library functions.  Only the bare minimum of IDL_ built-in
functions are implemented to support this.

.. _IDL: http://www.exelisvis.com/language/en-us/productsservices/idl.aspx
"""

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# ----------------------------------------------------------------------------

# For egg_info test builds to pass, put package imports here.
if not _ASTROPY_SETUP_:
    from .file_lines import file_lines
    from .median import median
    from .pcomp import pcomp
    from .rebin import rebin
    from .smooth import smooth
    from .uniq import uniq


class PydlException(Exception):
    """Base class for exceptions raised in PyDL functions.
    """
    pass

__all__ = ['file_lines', 'median', 'pcomp', 'rebin', 'smooth', 'uniq',
           'PydlException']
