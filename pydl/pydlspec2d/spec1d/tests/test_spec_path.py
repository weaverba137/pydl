# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-


def test_spec_path():
    import os
    from numpy import array
    from .. import spec_path
    if 'BOSS_SPECTRO_REDUX' in os.environ:
        bsr_orig = os.environ['BOSS_SPECTRO_REDUX']
        bsr = bsr_orig
    else:
        bsr_orig = None
        bsr = '/boss/spectro/redux'
        os.environ['BOSS_SPECTRO_REDUX'] = bsr
    if 'RUN2D' in os.environ:
        run2d_orig = os.environ['RUN2D']
        run2d = run2d_orig
    else:
        run2d_orig = None
        run2d = 'v1_2_3'
        os.environ['RUN2D'] = run2d
    p = spec_path(1234)
    assert p[0] == os.path.join(bsr, run2d, '1234')
    p = spec_path(1234, topdir=bsr, run2d=run2d)
    assert p[0] == os.path.join(bsr, run2d, '1234')
    p = spec_path(array([1234, 5678]), topdir=bsr, run2d=run2d)
    assert p[0] == os.path.join(bsr, run2d, '1234')
    assert p[1] == os.path.join(bsr, run2d, '5678')
    p = spec_path(1234, path=bsr)
    assert p[0] == bsr
    if bsr_orig is None:
        del os.environ['BOSS_SPECTRO_REDUX']
    if run2d_orig is None:
        del os.environ['RUN2D']
