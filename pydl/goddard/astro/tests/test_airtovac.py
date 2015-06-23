# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_airtovac():
    from .. import airtovac
    import numpy as np
    vacuum = airtovac(1900.0)
    assert vacuum == 1900.0
    vacuum = airtovac(2000.0)
    assert np.allclose(vacuum,2000.6475)
    air = np.array([1800.0,1900.0,2000.0,2100.0,2200.0,2300.0])
    vacuum = airtovac(air)
    assert np.allclose(vacuum,np.array([1800.0,1900.0,2000.6475,2100.6664,2200.6868,2300.7083]))
    vacuum = airtovac(6056.125)
    assert np.allclose(vacuum,6057.8019)
