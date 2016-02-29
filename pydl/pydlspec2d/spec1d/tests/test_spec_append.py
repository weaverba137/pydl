# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-


def test_spec_append():
    import os
    from numpy import array
    from .. import spec_append
    spec1 = array([[1, 1, 1, 1, 1],
                   [1, 1, 1, 1, 1]])
    spec2 = array([[2, 2, 2, 2, 2],
                   [2, 2, 2, 2, 2]])
    s = spec_append(spec1, spec2)
    assert (s == array([[1, 1, 1, 1, 1],
                        [1, 1, 1, 1, 1],
                        [2, 2, 2, 2, 2],
                        [2, 2, 2, 2, 2]])).all()
    spec2 = array([[2, 2, 2, 2],
                   [2, 2, 2, 2]])
    s = spec_append(spec1, spec2)
    assert (s == array([[1, 1, 1, 1, 1],
                        [1, 1, 1, 1, 1],
                        [2, 2, 2, 2, 0],
                        [2, 2, 2, 2, 0]])).all()
    s = spec_append(spec1, spec2, 1)
    assert (s == array([[1, 1, 1, 1, 1],
                        [1, 1, 1, 1, 1],
                        [0, 2, 2, 2, 2],
                        [0, 2, 2, 2, 2]])).all()
    spec1 = array([[1, 1, 1],
                   [1, 1, 1]])
    spec2 = array([[2, 2, 2, 2, 2],
                   [2, 2, 2, 2, 2]])
    s = spec_append(spec1, spec2, -2)
    assert (s == array([[0, 0, 1, 1, 1],
                        [0, 0, 1, 1, 1],
                        [2, 2, 2, 2, 2],
                        [2, 2, 2, 2, 2]])).all()
