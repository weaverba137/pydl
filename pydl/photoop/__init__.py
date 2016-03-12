# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
This subpackage implements functions from the photoop package.
"""
#
# Define this early on so that submodules can use it
#
from .. import PydlException


class PhotoopException(PydlException):
    """Exceptions raised by :mod:`pydl.photoop` that don't fit into a
    standard exception class like :exc:`ValueError`.
    """
    pass


__all__ = ['PhotoopException']
