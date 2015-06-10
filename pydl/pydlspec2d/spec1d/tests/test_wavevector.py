# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_wavevector():
    from .. import wavevector
    import numpy as np
    l = wavevector(3,5,wavemin=3,wavemax=5)
    ll = np.arange(20001,dtype='d') * 1.0e-04 + 3.0
    assert (l == ll).all()
