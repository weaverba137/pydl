# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-


def test_wavevector():
    from .. import wavevector
    import numpy as np
    l = wavevector(3, 4, binsz=0.1)
    ll = np.array([3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0])
    assert np.allclose(l, ll)
    l = wavevector(3, 4, wavemin=3, binsz=0.1)
    ll = np.array([3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0])
    assert np.allclose(l, ll)
