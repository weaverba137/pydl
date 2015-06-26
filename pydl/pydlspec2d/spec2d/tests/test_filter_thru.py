# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
#
def test_filter_thru():
    from .. import filter_thru
    import numpy as np
    from os.path import dirname, join
    from astropy.io import fits
    from astropy.tests.helper import raises
    # with fits.open(join('pydl','pydlspec2d','spec2d','tests','t','spPlate-4055-55359-0020.fits')) as hdulist:
    with fits.open(join(dirname(__file__),'t','spPlate-4055-55359-0020.fits')) as hdulist:
        flux = hdulist[0].data
        npix = hdulist[0].header['NAXIS1']
        ntrace = hdulist[0].header['NAXIS2']
        crval1 = hdulist[0].header['COEFF0']
        cd1_1 = hdulist[0].header['COEFF1']
    assert flux.shape == (ntrace,npix)
    loglam0 = crval1 + cd1_1*np.arange(npix,dtype=flux.dtype)
    waveimg = 10**(np.tile(loglam0,20).reshape(flux.shape))
    assert waveimg.shape == flux.shape
    f = filter_thru(flux,waveimg=waveimg)
    idl_data = np.array([
        [  14.1224,     0.127368,     0.508305,     0.293359,     0.381092,      2.51200,   -0.0306516,     0.298414,     0.185329,      1.11159,     0.166036,      1.08602,     0.162873,     0.529303,     0.290425,   -0.0399243,     0.209576,     0.778746,     0.578854,     0.429100,],
        [  18.3209,    0.0740831,     0.531655,     0.943787,     0.558309,      7.58337,     0.209387,     0.861768,     0.275096,      1.88644,     0.307924,      4.10726,     0.276175,     0.732941,     0.336644,    0.0456389,     0.725353,      2.25312,     0.413612,      1.08603,],
        [  14.4588,   -0.0113344,      1.15476,      3.76573,      1.35594,      17.3931,     0.647237,      2.96833,      1.00745,      1.46273,     0.815733,      11.2539,     0.957882,      2.57075,     0.821274,     0.162416,      2.08590,      5.34149,      1.32573,      3.67932,],
        [  10.4319,   -0.0250982,      1.83043,      5.57500,      1.64166,      19.9340,      1.18961,      3.89914,      1.83561,     0.884010,      1.22140,      13.7641,      1.66350,      3.35640,      1.41506,     0.105699,      3.04802,      6.47299,      2.31518,      5.10352,],
        [  7.59458,    0.0395837,      2.31145,      6.42070,      2.04900,      20.3303,      1.56850,      4.25935,      2.38152,     0.984143,      1.51820,      14.4626,      2.01851,      3.74417,      1.86966,     0.133837,      3.53355,      6.86765,      2.98000,      5.82073,]
        ],dtype=np.float32).T
    assert f.shape == (20,5)
    assert np.allclose(f,idl_data,atol=1.0e-7)
    #
    # Test bad input.
    #
    with raises(ValueError):
        f = filter_thru(flux)
    with raises(ValueError):
        f = filter_thru(flux,waveimg=waveimg,filter_prefix='sdss')
    return
#
#
#
def prepare_data():
    """Convert full spPlate file into a test-sized version.
    """
    from os import getenv
    from os.path import dirname, join
    from astropy.io import fits
    nTrace = 20
    spPlate = join(getenv('HOME'),'Downloads','spPlate-4055-55359.fits')
    spPlateOut =join(dirname(__file__),'t','spPlate-4055-55359-0020.fits')
    with fits.open(spPlate) as hdulist:
        newhdu = fits.PrimaryHDU(hdulist[0].data[0:20,:],header=hdulist[0].header)
        newhdulist = fits.HDUList([newhdu])
        newhdulist.writeto(spPlateOut)
    return 0
#
#
#
if __name__ == '__main__':
    from sys import exit
    exit(prepare_data())
