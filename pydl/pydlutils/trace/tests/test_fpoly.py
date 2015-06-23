# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_fpoly():
    from .. import fpoly
    import numpy as np
    from astropy.tests.helper import raises
    x = np.array([-1,-0.5,0,0.5,1],dtype='d')
    #
    # Test order
    #
    with raises(ValueError):
        f = fpoly(x,0)
    #
    # m = 1
    #
    f = fpoly(x,1)
    assert (f == np.ones((1,x.size),dtype='d')).all()
    #
    # m = 2
    #
    f = fpoly(x,2)
    foo = np.ones((2,x.size),dtype='d')
    foo[1,:] = x
    assert np.allclose(f,foo)
    #
    # m = 3
    #
    f = fpoly(x,3)
    foo = np.ones((3,x.size),dtype='d')
    foo[1,:] = x
    foo[2,:] = x**2
    assert np.allclose(f,foo)
    #
    # m = 4
    #
    f = fpoly(x,4)
    foo = np.ones((4,x.size),dtype='d')
    foo[1,:] = x
    foo[2,:] = x**2
    foo[3,:] = x**3
    assert np.allclose(f,foo)
    #
    # random float
    #
    f = fpoly(2.88,3)
    assert np.allclose(f,np.array([[1.00], [2.88], [8.29440]]))
