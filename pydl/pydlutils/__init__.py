# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
This subpackage implements functions from the idlutils package.
"""
#
# Define this early on so that submodules can use it
#
from astropy.utils.exceptions import AstropyUserWarning


class PydlutilsException(Exception):
    pass


class PydlutilsUserWarning(AstropyUserWarning):
    pass
