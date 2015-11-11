# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_cirrange():
    from .. import cirrange
    import numpy as np
    ra1 = np.linspace(-4.0*np.pi,4.0*np.pi,100)
    ra2 = cirrange(ra1,radians=True)
    assert (ra2 == (ra1 % (2.0*np.pi))).all()
    ra1 = np.rad2deg(ra1)
    ra2 = cirrange(ra1)
    assert (ra2 == (ra1 % 360.0)).all()
