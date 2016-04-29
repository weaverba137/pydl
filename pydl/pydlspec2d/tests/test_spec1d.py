# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
import os
# from astropy.tests.helper import raises
from ..spec1d import findspec, spec_append, spec_path, wavevector


class TestSpec1d(object):
    """Test the functions in pydl.pydlspec2d.spec1d.
    """

    def setup(self):
        # self.data_dir = join(dirname(__file__), 't')
        if 'BOSS_SPECTRO_REDUX' in os.environ:
            self.bsr_orig = os.environ['BOSS_SPECTRO_REDUX']
            self.bsr = bsr_orig
        else:
            self.bsr_orig = None
            self.bsr = '/boss/spectro/redux'
            os.environ['BOSS_SPECTRO_REDUX'] = self.bsr
        if 'RUN2D' in os.environ:
            self.run2d_orig = os.environ['RUN2D']
            self.run2d = self.run2d_orig
        else:
            self.run2d_orig = None
            self.run2d = 'v1_2_3'
            os.environ['RUN2D'] = self.run2d

    def teardown(self):
        if self.bsr_orig is None:
            del os.environ['BOSS_SPECTRO_REDUX']
        if self.run2d_orig is None:
            del os.environ['RUN2D']

    def test_findspec(self):
        """This is just a placeholder for now.
        """
        # slist = findspec(infile='file.in', sdss=True)
        assert True

    def test_spec_append(self):
        spec1 = np.array([[1, 1, 1, 1, 1],
                          [1, 1, 1, 1, 1]])
        spec2 = np.array([[2, 2, 2, 2, 2],
                          [2, 2, 2, 2, 2]])
        s = spec_append(spec1, spec2)
        assert (s == np.array([[1, 1, 1, 1, 1],
                               [1, 1, 1, 1, 1],
                               [2, 2, 2, 2, 2],
                               [2, 2, 2, 2, 2]])).all()
        spec2 = np.array([[2, 2, 2, 2],
                          [2, 2, 2, 2]])
        s = spec_append(spec1, spec2)
        assert (s == np.array([[1, 1, 1, 1, 1],
                               [1, 1, 1, 1, 1],
                               [2, 2, 2, 2, 0],
                               [2, 2, 2, 2, 0]])).all()
        s = spec_append(spec1, spec2, 1)
        assert (s == np.array([[1, 1, 1, 1, 1],
                               [1, 1, 1, 1, 1],
                               [0, 2, 2, 2, 2],
                               [0, 2, 2, 2, 2]])).all()
        spec1 = np.array([[1, 1, 1],
                          [1, 1, 1]])
        spec2 = np.array([[2, 2, 2, 2, 2],
                          [2, 2, 2, 2, 2]])
        s = spec_append(spec1, spec2, -2)
        assert (s == np.array([[0, 0, 1, 1, 1],
                               [0, 0, 1, 1, 1],
                               [2, 2, 2, 2, 2],
                               [2, 2, 2, 2, 2]])).all()

    def test_spec_path(self):
        p = spec_path(123)
        assert p[0] == os.path.join(self.bsr, self.run2d, '0123')
        p = spec_path(1234)
        assert p[0] == os.path.join(self.bsr, self.run2d, '1234')
        p = spec_path(1234, topdir=self.bsr, run2d=self.run2d)
        assert p[0] == os.path.join(self.bsr, self.run2d, '1234')
        p = spec_path(np.array([1234, 5678]), topdir=self.bsr,
                      run2d=self.run2d)
        assert p[0] == os.path.join(self.bsr, self.run2d, '1234')
        assert p[1] == os.path.join(self.bsr, self.run2d, '5678')
        p = spec_path(1234, path=self.bsr)
        assert p[0] == self.bsr

    def test_wavevector(self):
        l = wavevector(3, 4, binsz=0.1)
        ll = np.array([3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0])
        assert np.allclose(l, ll)
        l = wavevector(3, 4, wavemin=3, binsz=0.1)
        ll = np.array([3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0])
        assert np.allclose(l, ll)
