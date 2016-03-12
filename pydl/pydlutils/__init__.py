# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
This subpackage implements functions from the idlutils package.
"""
#
# Define this early on so that submodules can use it
#
from .. import PydlException
import astropy.utils.exceptions as aue


class PydlutilsException(PydlException):
    """Exceptions raised by :mod:`pydl.pydlutils` that don't fit into a
    standard exception class like :exc:`ValueError`.
    """
    pass


class PydlutilsUserWarning(aue.AstropyUserWarning):
    """Class for warnings issued by :mod:`pydl.pydlutils`.
    """
    pass


__all__ = ['PydlutilsException', 'PydlutilsUserWarning']
