# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import os
import numpy as np
from astropy.tests.helper import raises
from ..mangle import is_cap_used, read_fits_polygons


class TestMangle(object):
    """Test the functions in pydl.pydlutils.mangle.
    """

    def setup(self):
        self.data_dir = os.path.join(os.path.dirname(__file__), 't')

    def teardown(self):
        pass

    def test_is_cap_used(self):
        assert is_cap_used(1 << 2, 2)
        assert not is_cap_used(1 << 2, 1)

    def test_read_fits_polygons(self):
        poly = read_fits_polygons(os.path.join(self.data_dir, 'polygon.fits'))
        use_caps = np.array([31, 15, 31, 7, 31, 15, 15, 7, 15, 15,
                             15, 31, 15, 15, 15, 15, 15, 15, 31, 15],
                            dtype=np.uint32)
        assert (poly['USE_CAPS'] == use_caps).all()
        cm0 = np.array([-1.0, -0.99369437, 1.0, -1.0, 0.00961538])
        assert np.allclose(poly.CAPS.CM[0][0:poly.NCAPS[0]], cm0)
        assert poly[0]['NCAPS'] == 5


def fits_polygon_file():
    """Create a small test version of a FITS polygon file.
    """
    from datetime import date
    from sys import argv
    from astropy.io import fits
    from pydl import __version__ as pydlversion
    with fits.open(argv[1], uint=True) as hdulist:
        header0 = hdulist[0].header
        data = hdulist[1].data
    if 'DATE' in header0:
        header0['DATE'] = date.today().strftime('%Y-%m-%d')
    if 'IDLUTILS' in header0:
        header0['IDLUTILS'] = 'pydl-'+pydlversion
    hdu0 = fits.PrimaryHDU(header=header0)
    hdu1 = fits.BinTableHDU(data[0:20])
    hdulist2 = fits.HDUList([hdu0, hdu1])
    hdulist2.writeto('polygon.fits')
    return 0

if __name__ == '__main__':
    from sys import exit
    exit(fits_polygon_file())
