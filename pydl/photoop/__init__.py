# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
This subpackage implements functions from the photoop package.
"""
#
# Define this early on so that submodules can use it
#
#from .. import PydlException
#class PhotoopException(PydlException):
#    pass
class PhotoopException(Exception):
    pass
#
# Import subpackages
#
from . import photoobj
from . import sdssio
from . import window
#
# Set __all__
#
__all__ = ['PhotoopException', 'photoobj', 'sdssio', 'window' ]
