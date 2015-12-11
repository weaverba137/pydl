# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
This subpackage implements functions from the idlutils package.
"""
#
# Define this early on so that submodules can use it
#
import astropy.utils.exceptions as aue


class PydlutilsException(Exception):
    pass


class PydlutilsUserWarning(aue.AstropyUserWarning):
    pass
