# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import os
import numpy as np
from astropy.tests.helper import raises
from .. import PydlutilsException
from ..mangle import (ManglePolygon, is_cap_used, read_fits_polygons,
                      read_mangle_polygons, set_use_caps)


class TestMangle(object):
    """Test the functions in pydl.pydlutils.mangle.
    """

    def setup(self):
        self.data_dir = os.path.join(os.path.dirname(__file__), 't')
        self.poly_fits = os.path.join(self.data_dir, 'polygon.fits')
        self.poly_ply = os.path.join(self.data_dir, 'polygon.ply')
        self.bad_ply = os.path.join(self.data_dir, 'median_data.txt')

    def teardown(self):
        pass

    def test_ManglePolygon(self):
        with raises(ValueError):
            poly = ManglePolygon(weight=1.0)
        with raises(ValueError):
            poly = ManglePolygon()
        x = np.array([[0.0, 0.0, 1.0],
                      [1.0, 0.0, 0.0],
                      [0.0, 1.0, 1.0]])
        cm = np.array([1.0, 1.0, 1.0])
        poly = ManglePolygon(x=x, cm=cm, str=np.pi/2.0)
        assert poly.ncaps == 3
        assert poly.weight == 1.0
        assert poly.use_caps == (1 << 3) - 1
        poly = ManglePolygon(x=x, cm=cm, weight=0.5)
        assert poly.weight == 0.5
        poly = ManglePolygon(x=x, cm=cm, pixel=20)
        assert poly.pixel == 20
        poly = ManglePolygon(x=x, cm=cm, use_caps=3, str=np.pi/2.0)
        assert poly.use_caps == 3
        poly2 = poly.copy()
        assert poly2.use_caps == poly.use_caps
        assert (poly2.cm == poly.cm).all()
        assert np.allclose(poly2.str, np.pi/2.0)

    def test_is_cap_used(self):
        assert is_cap_used(1 << 2, 2)
        assert not is_cap_used(1 << 2, 1)

    def test_read_fits_polygons(self):
        poly = read_fits_polygons(self.poly_fits)
        use_caps = np.array([31, 15, 31, 7, 31, 15, 15, 7, 15, 15,
                             15, 31, 15, 15, 15, 15, 15, 15, 31, 15],
                            dtype=np.uint32)
        #
        # Attribute access doesn't work on unsigned columns.
        #
        assert (poly['USE_CAPS'] == use_caps).all()
        assert (poly['use_caps'] == use_caps).all()
        with raises(AttributeError):
            foo = poly.no_such_attribute
        cm0 = np.array([-1.0, -0.99369437, 1.0, -1.0, 0.00961538])
        assert np.allclose(poly.cm[0][0:poly.ncaps[0]], cm0)
        assert poly[0]['NCAPS'] == 5
        poly = read_fits_polygons(self.poly_fits, convert=True)
        assert poly[0].use_caps == 31
        assert np.allclose(poly[0].cm, cm0)
        assert poly[0].cmminf == 4

    def test_read_mangle_polygons(self):
        with raises(PydlutilsException):
            poly = read_mangle_polygons(self.bad_ply)
        poly = read_mangle_polygons(self.poly_ply)
        assert len(poly.header) == 3
        assert poly.header[0] == 'pixelization 6s'
        assert len(poly) == 4
        assert np.allclose(poly[0].x[0, :],
                           np.array([0.0436193873653360, 0.9990482215818578,
                                     0.0]))
        assert poly[3].ncaps == 3

    def test_set_use_caps(self):
        poly = read_fits_polygons(self.poly_fits, convert=True)
        old_use_caps = poly[0].use_caps
        index_list = list(range(poly[0].ncaps))
        use_caps = set_use_caps(poly[0], index_list, allow_doubles=True)
        assert use_caps == poly[0].use_caps
        use_caps = set_use_caps(poly[0], index_list)
        assert use_caps == poly[0].use_caps


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
