# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
This subpackage corresponds to the sdss directory of idlutils.
"""
from sdss_flagexist import sdss_flagexist
from sdss_flagname import sdss_flagname
from sdss_flagval import sdss_flagval
from set_maskbits import set_maskbits
#
# Cache the maskbits file.
#
maskbits = set_maskbits()
