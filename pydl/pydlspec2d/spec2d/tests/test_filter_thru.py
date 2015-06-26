# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
#
def test_filter_thru():
    assert True
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
