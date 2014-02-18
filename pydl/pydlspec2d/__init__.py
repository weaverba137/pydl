# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
This subpackage implements functions from the photoop package.
"""
#
# Define this early on so that submodules can use it
#
#from .. import PydlException
#class Pydlspec2dException(PydlException):
#    pass
class Pydlspec2dException(Exception):
    pass
#
# Import subpackages
#
from . import spec1d
from . import spec2d
#
# Set __all__
#
__all__ = ['Pydlspec2dException', 'spec1d', 'spec2d']
