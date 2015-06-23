# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_vactoair():
    from .. import vactoair
    import numpy as np
    air = vactoair(1900.0)
    assert air == 1900.0
    air = vactoair(2000.0)
    assert np.allclose(air,1999.3526)
    vacuum = np.array([1800.0,1900.0,2000.0,2100.0,2200.0,2300.0])
    air = vactoair(vacuum)
    assert np.allclose(air,np.array([1800.0,1900.0,1999.3526,2099.3337,2199.3133,2299.2918]))
