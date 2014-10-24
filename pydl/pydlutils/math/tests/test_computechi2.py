# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_computechi2():
    import numpy as np
    from .. import computechi2
    x = np.arange(20)
    y = np.array([6.6198438, 1.3491303, 0.41035045, 9.4013375, 4.1103360,
        4.3522868, 4.6338078, 4.7400367, 5.1860726, 5.1082748,
        5.1084596, 5.2990997, 5.5987537, 5.7007504, 5.7855296,
        4.1123709, 6.9437957, 4.9956179, 4.3724215, 3.6245063])
    sqivar = np.array([1.32531, 0.886090, 1.08384, 1.04489, 1.46807,
        1.30800, 0.507725, 1.12840, 0.955025, 1.35925,
        1.10126, 1.45690, 0.575700, 0.949710, 1.23368,
        0.536489, 0.772543, 0.957729, 0.883976, 1.11559])
    templates = np.vstack((np.ones((20,),dtype='d'),x)).transpose()
    chi2 = computechi2(y,sqivar,templates)
    #
    # 20 data points, 2 parameters = 18 degrees of freedom.
    #
    assert chi2.dof == 18
    #
    # The covariance matrix should be symmetric
    #
    assert (chi2.covar.T == chi2.covar).all()
    #
    # The variance vector should be the same as the diagonal of the covariance
    # matrix
    #
    assert (chi2.var == np.diag(chi2.covar)).all()
