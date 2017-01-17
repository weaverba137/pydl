# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
import os
from astropy.tests.helper import raises
from astropy.utils.data import get_pkg_data_filename
from .. import Pydlspec2dException
from ..spec1d import (HMF, findspec, spec_append, spec_path, template_metadata,
                      wavevector)


class TestSpec1d(object):
    """Test the functions in pydl.pydlspec2d.spec1d.
    """

    def setup(self):
        self.env = {'BOSS_SPECTRO_REDUX': '/boss/spectro/redux',
                    'SPECTRO_REDUX': '/sdss/spectro/redux',
                    'RUN2D': 'v1_2_3',
                    'RUN1D': 'v1_2_3'}
        self.original_env = dict()
        for key in self.env:
            if key in os.environ:
                self.original_env[key] = os.environ[key]
            else:
                self.original_env[key] = None
            os.environ[key] = self.env[key]

    def teardown(self):
        for key in self.original_env:
            if self.original_env[key] is None:
                del os.environ[key]
            else:
                os.environ[key] = self.original_env[key]

    def test_findspec(self):
        """This is just a placeholder for now.
        """
        # slist = findspec(infile='file.in', sdss=True)
        assert True

    def test_hmf_init(self):
        """Test initialization of HMF object
        """
        spec = np.random.random((20, 100))
        invvar = np.random.random((20, 100))
        hmf = HMF(spec, invvar)
        assert hmf.K == 4
        assert hmf.log.level == 20  # INFO
        hmf = HMF(spec, invvar, K=6, verbose=True)
        assert hmf.K == 6
        assert hmf.log.level == 10  # DEBUG

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
        bsr = self.env['BOSS_SPECTRO_REDUX']
        run2d = self.env['RUN2D']
        p = spec_path(123)
        assert p[0] == os.path.join(bsr, run2d, '0123')
        p = spec_path(1234)
        assert p[0] == os.path.join(bsr, run2d, '1234')
        p = spec_path(1234, topdir=bsr, run2d=run2d)
        assert p[0] == os.path.join(bsr, run2d, '1234')
        p = spec_path(np.array([1234, 5678]), topdir=bsr, run2d=run2d)
        assert p[0] == os.path.join(bsr, run2d, '1234')
        assert p[1] == os.path.join(bsr, run2d, '5678')
        p = spec_path(1234, path=bsr)
        assert p[0] == bsr

    def test_template_metadata(self):
        with raises(Pydlspec2dException):
            slist, metadata = template_metadata('/no/such/file.par')
        inputfile = get_pkg_data_filename('t/test_template_metadata.par')
        slist, metadata = template_metadata(inputfile)
        assert metadata['object'] == 'gal'
        assert not metadata['nonnegative']

    def test_wavevector(self):
        l = wavevector(3, 4, binsz=0.1)
        ll = np.array([3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0])
        assert np.allclose(l, ll)
        l = wavevector(3, 4, wavemin=3, binsz=0.1)
        ll = np.array([3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0])
        assert np.allclose(l, ll)
