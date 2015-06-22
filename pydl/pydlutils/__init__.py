# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
This subpackage implements functions from the idlutils package.
"""
#
# Define this early on so that submodules can use it
#
#from .. import PydlException
#class PydlutilsException(PydlException):
#    pass
class PydlutilsException(Exception):
    pass
#
from astropy.utils.exceptions import AstropyUserWarning
class PydlutilsUserWarning(AstropyUserWarning):
    pass
#
# Import subpackages
#
from . import bspline
from . import cooling
#from . import coord
from . import image
from . import mangle
from . import math
from . import misc
from . import sdss
from . import spheregroup
from . import trace
from . import yanny
#
# Set __all__
#
__all__ = ['PydlutilsException', 'bspline', 'cooling', 'image', 'mangle',
    'math', 'misc', 'sdss', 'spheregroup', 'trace', 'yanny']
