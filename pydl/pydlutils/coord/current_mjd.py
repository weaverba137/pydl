# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def current_mjd():
    """Return the current MJD.
    """
    from ...goddard.astro import get_juldate
    return get_juldate() - 2400000.5
