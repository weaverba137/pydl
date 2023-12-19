# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
====
pydl
====

Python replacements for functions that are part of the `IDL®`_ built-in library, or
part of astronomical `IDL®`_ libraries.  The emphasis is on reproducing results of
the astronomical library functions.  Only the bare minimum of `IDL®`_ built-in
functions are implemented to support this.

.. _`IDL®`: https://www.nv5geospatialsoftware.com/Products/IDL
"""
import os
from astropy.tests.runner import TestRunner

try:
    from .version import version as __version__
except ImportError:
    __version__ = ''

from .file_lines import file_lines
from .median import median
from .pcomp import pcomp
from .rebin import rebin
from .smooth import smooth
from .uniq import uniq


# Create the test function for self test
test = TestRunner.make_test_runner_in(os.path.dirname(__file__))
test.__test__ = False


class PydlException(Exception):
    """Base class for exceptions raised in PyDL functions.
    """
    pass


__all__ = ['file_lines', 'median', 'pcomp', 'rebin', 'smooth', 'uniq',
           'PydlException']

# Clean up namespace
del os
del TestRunner
