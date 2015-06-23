# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_fchebyshev_split():
    from .. import fchebyshev_split
    import numpy as np
    from astropy.tests.helper import raises
    x = np.array([-1,-0.5,0,0.5,1],dtype='d')
    #
    # Test order
    #
    with raises(ValueError):
        f = fchebyshev_split(x,0)
    with raises(ValueError):
        f = fchebyshev_split(x,1)
    #
    # m = 2
    #
    f = fchebyshev_split(x,2)
    foo = np.ones((2,x.size),dtype='d')
    foo[0,:] = (x >= 0).astype(x.dtype)
    assert np.allclose(f,foo)
    #
    # m = 3
    #
    f = fchebyshev_split(x,3)
    foo = np.ones((3,x.size),dtype='d')
    foo[0,:] = (x >= 0).astype(x.dtype)
    foo[2,:] = x
    assert np.allclose(f,foo)
    #
    # m = 4
    #
    f = fchebyshev_split(x,4)
    foo = np.ones((4,x.size),dtype='d')
    foo[0,:] = (x >= 0).astype(x.dtype)
    foo[2,:] = x
    foo[3,:] = (2.0*x**2 - 1.0)
    assert np.allclose(f,foo)
    #
    # m = 5
    #
    f = fchebyshev_split(x,5)
    foo = np.ones((5,x.size),dtype='d')
    foo[0,:] = (x >= 0).astype(x.dtype)
    foo[2,:] = x
    foo[3,:] = (2.0*x**2 - 1.0)
    foo[4,:] = (4.0*x**3 - 3.0*x)
    assert np.allclose(f,foo)
    #
    # random float
    #
    f = fchebyshev_split(2.88,3)
    assert np.allclose(f,np.array([[1.00], [1.00], [2.88]]))
