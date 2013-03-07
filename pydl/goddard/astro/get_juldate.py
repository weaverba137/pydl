# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from __future__ import print_function
#
#
#
def get_juldate():
    """Returns the current Julian date.

    Uses the MJD trick & adds the offset to get JD.
    """
    import time
    mjd = time.time()/86400.0 + 40587.0
    return mjd + 2400000.5
#
#
#
def main(args=None):
    """Allow this module to be run in scripts.
    """
    jd = get_juldate()
    print(jd)
    return
