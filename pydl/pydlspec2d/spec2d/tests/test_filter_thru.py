# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-


def test_filter_thru():
    from .. import filter_thru
    import numpy as np
    from os.path import dirname, join
    from astropy.io import fits
    from astropy.tests.helper import raises
    # with fits.open(join('pydl','pydlspec2d','spec2d','tests','t','spPlate-4055-55359-0020.fits')) as hdulist:
    with fits.open(join(dirname(__file__), 't', 'spPlate-4055-55359-0020.fits')) as hdulist:
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
    idl_data_file = join(dirname(__file__), 't', 'filter_thru_idl_data.txt')
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
    from os import getenv
    from os.path import dirname, join
    from astropy.io import fits
    nTrace = 20
    spPlate = join(getenv('HOME'), 'Downloads', 'spPlate-4055-55359.fits')
    spPlateOut = join(dirname(__file__), 't', 'spPlate-4055-55359-0020.fits')
    with fits.open(spPlate) as hdulist:
        newhdu = fits.PrimaryHDU(hdulist[0].data[0:20, :],
                                 header=hdulist[0].header)
        newhdulist = fits.HDUList([newhdu])
        newhdulist.writeto(spPlateOut)
    return 0


if __name__ == '__main__':
    from sys import exit
    exit(prepare_data())
