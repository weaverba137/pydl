# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
from astropy.tests.helper import raises
from astropy.utils.data import get_pkg_data_filename
from .. import PydlutilsException
from .. import mangle as mng


class TestMangle(object):
    """Test the functions in pydl.pydlutils.mangle.
    """

    def setup(self):
        self.poly_fits = get_pkg_data_filename('t/polygon.fits')
        self.no_id_fits = get_pkg_data_filename('t/polygon_no_id.fits')
        self.one_cap_fits = get_pkg_data_filename('t/polygon_one_cap.fits')
        self.poly_ply = get_pkg_data_filename('t/polygon.ply')
        self.bad_ply = get_pkg_data_filename('t/median_data.txt')

    def teardown(self):
        pass

    def test_ManglePolygon(self):
        #
        # Zero caps
        #
        poly = mng.ManglePolygon()
        assert np.allclose(poly.str, 4.0*np.pi)
        assert poly.cmminf() is None
        assert not poly.gzeroar()
        assert np.allclose(poly.garea(), 4.0*np.pi)
        #
        # One cap.
        #
        x = np.array([[0.0, 0.0, 1.0]])
        cm = np.array([0.5])
        poly = mng.ManglePolygon(x=x, cm=cm)
        assert np.allclose(poly.garea(), np.pi)
        #
        # Bad inputs
        #
        with raises(ValueError):
            poly = mng.ManglePolygon(weight=1.0)
        #
        # Multiple caps
        #
        x = np.array([[0.0, 0.0, 1.0],
                      [1.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0]])
        cm = np.array([1.0, 1.0, 1.0])
        poly = mng.ManglePolygon(x=x, cm=cm, str=np.pi/2.0)
        assert poly.ncaps == 3
        assert poly.weight == 1.0
        assert poly.use_caps == (1 << 3) - 1
        poly = mng.ManglePolygon(x=x, cm=cm, weight=0.5)
        assert poly.weight == 0.5
        poly = mng.ManglePolygon(x=x, cm=cm, pixel=20)
        assert poly.pixel == 20
        poly = mng.ManglePolygon(x=x, cm=cm, use_caps=3, str=np.pi/2.0)
        assert poly.use_caps == 3
        poly2 = poly.copy()
        assert poly2.use_caps == poly.use_caps
        assert (poly2.cm == poly.cm).all()
        assert np.allclose(poly2.str, np.pi/2.0)
        x = np.array([[0.0, 0.0, 1.0],
                      [1.0, 0.0, 0.0]])
        cm = np.array([1.0, 1.0])
        poly = mng.ManglePolygon(x=x, cm=cm)
        poly2 = poly.add_caps(np.array([[0.0, 1.0, 0.0], ]), np.array([1.0, ]))
        assert poly2.ncaps == 3
        assert poly2.use_caps == poly.use_caps
        assert poly2.str == 1.0  # dummy value!
        poly3 = poly.polyn(poly2, 2)
        assert poly3.ncaps == 3
        assert poly3.use_caps == poly.use_caps
        assert np.allclose(poly3.x[2, :], np.array([0.0, 1.0, 0.0]))
        poly3 = poly.polyn(poly2, 2, complement=True)
        assert poly3.ncaps == 3
        assert poly3.use_caps == poly.use_caps
        assert np.allclose(poly3.cm[2], -1.0)

    def test_angles_to_x(self):
        x = mng.angles_to_x(np.array([[0.0, 0.0], [90.0, 90.0],
                                      [0.0, 90.0]]))
        assert np.allclose(x, np.array([[0.0, 0.0, 1.0], [0.0, 1.0, 0.0],
                                        [1.0, 0.0, 0.0]]))
        a = mng.angles_to_x(np.array([[0.0, 90.0], [90.0, 0.0],
                                      [0.0, 0.0]]), latitude=True)
        assert np.allclose(x, np.array([[0.0, 0.0, 1.0], [0.0, 1.0, 0.0],
                                        [1.0, 0.0, 0.0]]))

    def test_cap_distance(self):
        x = np.array([0.0, 0.0, 1.0])
        cm = 1.0
        with raises(ValueError):
            d = mng.cap_distance(x, cm, np.array([[1.0, 2.0, 3.0, 4.0], ]))
        d = mng.cap_distance(x, cm, np.array([[0.0, 45.0], ]))
        assert np.allclose(d, np.array([45.0]))
        y = mng.angles_to_x(np.array([[0.0, 45.0], ]), latitude=True)
        d = mng.cap_distance(x, cm, y)
        assert np.allclose(d, np.array([45.0]))
        d = mng.cap_distance(x, cm, np.array([[0.0, -45.0], ]))
        assert np.allclose(d, np.array([-45.0]))
        d = mng.cap_distance(x, -1.0, np.array([[0.0, -45.0], ]))
        assert np.allclose(d, np.array([45.0]))

    def test_circle_cap(self):
        with raises(ValueError):
            x, cm = mng.circle_cap(90.0, np.array([[1.0, 2.0, 3.0, 4.0], ]))
        xin = np.array([[0.0, 0.0, 1.0], ])
        x, cm = mng.circle_cap(90.0, xin)
        assert np.allclose(x, xin)
        assert np.allclose(cm, 1.0)
        radec = np.array([[0.0, 90.0], ])
        x, cm = mng.circle_cap(90.0, radec)
        assert np.allclose(x, xin)
        assert np.allclose(cm, 1.0)
        x, cm = mng.circle_cap(np.float32(90.0), radec)
        assert np.allclose(x, xin)
        assert np.allclose(cm, 1.0)
        x, cm = mng.circle_cap(np.array([90.0, ]), radec)
        assert np.allclose(x, xin)
        assert np.allclose(cm, np.array([1.0, ]))
        with raises(ValueError):
            x, cm = mng.circle_cap(np.array([90.0, 90.0]), radec)

    def test_is_cap_used(self):
        assert mng.is_cap_used(1 << 2, 2)
        assert not mng.is_cap_used(1 << 2, 1)

    def test_is_in_cap(self):
        x = np.array([0.0, 0.0, 1.0])
        cm = 1.0
        d = mng.is_in_cap(x, cm, np.array([[0.0, 45.0], [0.0, -45.0]]))
        assert (d == np.array([True, False])).all()

    def test_is_in_polygon(self):
        x = np.array([[0.0, 0.0, 1.0],
                      [0.0, 1.0, 0.0],
                      [1.0, 0.0, 0.0]])
        cm = np.array([1.0, 1.0, 1.0])
        p = mng.ManglePolygon(x=x, cm=cm)
        d = mng.is_in_polygon(p, np.array([[45.0, 45.0], [135.0, -45.0]]))
        assert (d == np.array([True, False])).all()
        d = mng.is_in_polygon(p, np.array([[45.0, 45.0], [-45.0, 45.0]]),
                              ncaps=2)
        assert (d == np.array([True, False])).all()
        poly = mng.read_fits_polygons(self.no_id_fits)
        d = mng.is_in_polygon(poly[0],
                              np.array([[0.0, 10.0],
                                        [4.0, 4.5],
                                        [90, 4.0],
                                        [135, 3.0],
                                        [180, 0.5],
                                        [270, 0.25]]))
        assert (d == np.array([False, True, False, False, False, False])).all()

    def test_is_in_window(self):
        np.random.seed(271828)
        RA = 7.0*np.random.random(1000) + 268.0
        Dec = 90.0 - np.degrees(np.arccos(0.08*np.random.random(1000)))
        points = np.vstack((RA, Dec)).T
        poly = mng.read_fits_polygons(self.poly_fits)
        i = mng.is_in_window(poly, points)
        assert i[0].sum() == 3

    def test_read_fits_polygons(self):
        poly = mng.read_fits_polygons(self.poly_fits)
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
        poly = mng.read_fits_polygons(self.poly_fits, convert=True)
        assert poly[0].use_caps == 31
        assert np.allclose(poly[0].cm, cm0)
        assert poly[0].cmminf() == 4
        #
        # A FITS file might not contain IFIELD.
        #
        poly = mng.read_fits_polygons(self.no_id_fits)
        assert len(poly) == 1
        #
        # A FITS file might contain exactly one polygon with exactly one cap.
        #
        poly = mng.read_fits_polygons(self.one_cap_fits, convert=True)
        assert poly[0].ncaps == 1

    def test_read_mangle_polygons(self):
        with raises(PydlutilsException):
            poly = mng.read_mangle_polygons(self.bad_ply)
        poly = mng.read_mangle_polygons(self.poly_ply)
        assert len(poly.header) == 3
        assert poly.header[0] == 'pixelization 6s'
        assert len(poly) == 4
        assert np.allclose(poly[0].x[0, :],
                           np.array([0.0436193873653360, 0.9990482215818578,
                                     0.0]))
        assert poly[3].ncaps == 3

    def test_set_use_caps(self):
        poly = mng.read_fits_polygons(self.poly_fits, convert=True)
        old_use_caps = poly[0].use_caps
        index_list = list(range(poly[0].ncaps))
        use_caps = mng.set_use_caps(poly[0], index_list, allow_doubles=True)
        assert use_caps == poly[0].use_caps
        use_caps = mng.set_use_caps(poly[0], index_list)
        assert use_caps == poly[0].use_caps
        x = np.array([[0.0, 0.0, 1.0],
                      [0.0, 1.0, 0.0],
                      [1.0, 0.0, 0.0],
                      [1.0-1.0e-8, 1.0e-8, 0.0]])
        cm = np.array([1.0, 1.0, 1.0, 1.0-1.0e-8])
        p = mng.ManglePolygon(x=x, cm=cm)
        index_list = list(range(p.ncaps))
        assert p.use_caps == 2**4 - 1
        use_caps = mng.set_use_caps(p, index_list, tol=1.0e-7)
        assert use_caps == 2**3 - 1
        x = np.array([[0.0, 0.0, 1.0],
                      [0.0, 1.0, 0.0],
                      [1.0, 0.0, 0.0],
                      [1.0-1.0e-8, 1.0e-8, 0.0]])
        cm = np.array([1.0, 1.0, 1.0, -1.0+1.0e-8])
        p = mng.ManglePolygon(x=x, cm=cm)
        index_list = list(range(p.ncaps))
        assert p.use_caps == 2**4 - 1
        use_caps = mng.set_use_caps(p, index_list, tol=1.0e-7)
        assert use_caps == 2**3 - 1
        use_caps = mng.set_use_caps(p, index_list, tol=1.0e-7,
                                    allow_neg_doubles=True)
        assert use_caps == 2**4 - 1

    def test_x_to_angles(self):
        a = mng.x_to_angles(np.array([[0.0, 0.0, 1.0], [0.0, 1.0, 0.0],
                                      [1.0, 0.0, 0.0]]))
        assert np.allclose(a, np.array([[0.0, 0.0], [90.0, 90.0],
                                        [0.0, 90.0]]))
        a = mng.x_to_angles(np.array([[0.0, 0.0, 1.0], [0.0, 1.0, 0.0],
                                      [1.0, 0.0, 0.0]]), latitude=True)
        assert np.allclose(a, np.array([[0.0, 90.0], [90.0, 0.0],
                                        [0.0, 0.0]]))

    def test_single_polygon(self):
        x = np.array([[0.0, 0.0, 1.0],
                      [0.0, 1.0, 0.0],
                      [1.0, 0.0, 0.0],
                      [1.0-1.0e-8, 1.0e-8, 0.0]])
        cm = np.array([1.0, 1.0, 1.0, -1.0+1.0e-8])
        p = mng.ManglePolygon(x=x, cm=cm)
        p2 = mng._single_polygon(p)
        assert id(p2) == id(p)
        pl = mng.PolygonList()
        pl.append(p)
        p2 = mng._single_polygon(pl)
        assert id(p2) == id(p)
        with raises(ValueError):
            pl.append(p2)
            p3 = mng._single_polygon(pl)
        p = mng.read_fits_polygons(self.no_id_fits)
        p2 = mng._single_polygon(p)
        assert p2.ncaps == 4
        p2 = mng._single_polygon(p[0])
        assert p2.use_caps == 15


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
