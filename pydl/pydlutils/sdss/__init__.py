# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
This subpackage corresponds to the sdss directory of idlutils.
"""
from .default_skyversion import default_skyversion
from .sdss_astrombad import sdss_astrombad
from .sdss_flagexist import sdss_flagexist
from .sdss_flagname import sdss_flagname
from .sdss_flagval import sdss_flagval
from .sdss_objid import sdss_objid
from .sdss_sweep_circle import sdss_sweep_circle
#
# Cache for the maskbits file.
#
maskbits = None
#
# Cache the sweep index
#
sweep_cache = {'star':None,'gal':None,'sky':None}
#
# Cache sdss_astrombad data
#
opbadfields = None
