# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_cirrange():
    from .. import cirrange
    import numpy as np
    ra1 = np.linspace(-4.0*np.pi,4.0*np.pi,100)
    ra2 = cirrange(ra2,radians=True)
    assert ra2 == ra1 % 2.0*np.pi
