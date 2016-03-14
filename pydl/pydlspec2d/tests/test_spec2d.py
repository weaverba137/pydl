# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
import os
from astropy.io import fits
from astropy.tests.helper import raises
from ..spec2d import aesthetics, combine1fiber, filter_thru
from .. import Pydlspec2dException


class TestSpec2d(object):
    """Test the functions in pydl.pydlspec2d.spec2d.
    """

    def setup(self):
        self.data_dir = os.path.join(os.path.dirname(__file__), 't')
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

    def test_aesthetics(self):
        np.random.seed(137)
        flux = np.random.rand(100)
        ivar = np.random.rand(100)
        #
        # No bad
        #
        f = aesthetics(flux, ivar)
        assert (f == flux).all()
        #
        # Bad points
        #
        ivar[ivar < 0.1] = 0.0
        #
        # Bad method
        #
        with raises(Pydlspec2dException):
            f = aesthetics(flux, ivar, 'badmethod')
        #
        # Nothing
        #
        f = aesthetics(flux, ivar, 'nothing')
        assert (f == flux).all()

    def test_combine1fiber(self):
        pass

    def test_filter_thru(self):
        with fits.open(os.path.join(self.data_dir,
                                    'spPlate-4055-55359-0020.fits')) as hdulist:
            flux = hdulist[0].data
            npix = hdulist[0].header['NAXIS1']
            ntrace = hdulist[0].header['NAXIS2']
            crval1 = hdulist[0].header['COEFF0']
            cd1_1 = hdulist[0].header['COEFF1']
        assert flux.shape == (ntrace, npix)
        loglam0 = crval1 + cd1_1*np.arange(npix, dtype=flux.dtype)
        waveimg = 10**(np.tile(loglam0, 20).reshape(flux.shape))
        assert waveimg.shape == flux.shape
        f = filter_thru(flux, waveimg=waveimg)
        idl_data_file = os.path.join(self.data_dir,
                                     'filter_thru_idl_data.txt')
        idl_data = np.loadtxt(idl_data_file, dtype='f', delimiter=',').T
        assert f.shape == (20, 5)
        assert np.allclose(f, idl_data, atol=1.0e-6)
        #
        # Test bad input.
        #
        with raises(ValueError):
            f = filter_thru(flux)
        with raises(ValueError):
            f = filter_thru(flux, waveimg=waveimg, filter_prefix='sdss')
        return


def prepare_data():
    """Convert full spPlate file into a test-sized version.
    """
    nTrace = 20
    spPlate = os.path.join(os.getenv('HOME'), 'Downloads',
                           'spPlate-4055-55359.fits')
    spPlateOut = os.path.join(os.path.dirname(__file__), 't',
                              'spPlate-4055-55359-0020.fits')
    with fits.open(spPlate) as hdulist:
        newhdu = fits.PrimaryHDU(hdulist[0].data[0:20, :],
                                 header=hdulist[0].header)
        newhdulist = fits.HDUList([newhdu])
        newhdulist.writeto(spPlateOut)
    return 0


if __name__ == '__main__':
    from sys import exit
    exit(prepare_data())
