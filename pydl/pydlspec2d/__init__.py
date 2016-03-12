# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
This subpackage implements functions from the idlspec2d package.
"""
#
# Define this early on so that submodules can use it
#
from .. import PydlException
import astropy.utils.exceptions as aue


class Pydlspec2dException(PydlException):
    """Exceptions raised by :mod:`pydl.pydlspec2d` that don't fit into a
    standard exception class like :exc:`ValueError`.
    """
    pass


class Pydlspec2dUserWarning(aue.AstropyUserWarning):
    """Class for warnings issued by :mod:`pydl.pydlspec2d`.
    """
    pass


__all__ = ['Pydlspec2dException', 'Pydlspec2dUserWarning']
