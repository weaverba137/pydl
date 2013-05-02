# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
This subpackage corresponds to the sdss directory of idlutils.
"""
from default_skyversion import default_skyversion
from sdss_flagexist import sdss_flagexist
from sdss_flagname import sdss_flagname
from sdss_flagval import sdss_flagval
from sdss_objid import sdss_objid
from sdss_sweep_circle import sdss_sweep_circle
from set_maskbits import set_maskbits
#
# Cache the maskbits file.
#
maskbits = set_maskbits()
#
# Remove this from the namespace after use.
#
del set_maskbits
#
# Cache the sweep index
#
sweep_cache = {'star':None,'gal':None,'sky':None}
