# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from os import remove
import numpy as np
import tempfile
from astropy.tests.helper import raises
from .. import PydlutilsException
from ..misc import djs_laxisgen, djs_laxisnum, hogg_iau_name, struct_print


class TestMisc(object):
    """Test the functions in pydl.pydlutils.misc.
    """

    def setup(self):
        pass

    def teardown(self):
        pass

    def test_djs_laxisgen(self):
        #
        # 1d
        #
        assert (np.arange(4, dtype='i4') == djs_laxisgen((4,))).all()
        #
        # 2d
        #
        l = np.array([[0, 0, 0, 0],
                        [1, 1, 1, 1],
                        [2, 2, 2, 2],
                        [3, 3, 3, 3]],
                        dtype='i4')
        assert (l == djs_laxisgen((4, 4))).all()
        assert (l.T == djs_laxisgen((4, 4), iaxis=1)).all()
        with raises(ValueError):
            foo = djs_laxisgen((4, 4), iaxis=2)
        #
        # 3d
        #
        l = np.zeros((4, 4, 4), dtype='i4')
        l[1, :, :] = 1
        l[2, :, :] = 2
        l[3, :, :] = 3
        assert (l == djs_laxisgen((4, 4, 4))).all()
        assert (l.swapaxes(0, 1) == djs_laxisgen((4, 4, 4), iaxis=1)).all()
        assert (l.swapaxes(0, 2) == djs_laxisgen((4, 4, 4), iaxis=2)).all()
        with raises(ValueError):
            foo = djs_laxisgen((4, 4, 4), iaxis=3)
        #
        # More d
        #
        with raises(ValueError):
            foo = djs_laxisgen((4, 4, 4, 4))

    def test_djs_laxisnum(self):
        #
        # 1d
        #
        assert (np.zeros((4,), dtype='i4') == djs_laxisnum((4,))).all()
        #
        # 2d
        #
        l = np.array([[0, 0, 0, 0],
                        [1, 1, 1, 1],
                        [2, 2, 2, 2],
                        [3, 3, 3, 3]],
                        dtype='i4')
        assert (l == djs_laxisnum((4, 4))).all()
        assert (l.T == djs_laxisnum((4, 4), iaxis=1)).all()
        with raises(ValueError):
            foo = djs_laxisnum((4, 4), iaxis=2)
        #
        # 3d
        #
        l = np.zeros((4, 4, 4), dtype='i4')
        l[1, :, :] = 1
        l[2, :, :] = 2
        l[3, :, :] = 3
        assert (l == djs_laxisnum((4, 4, 4))).all()
        assert (l.swapaxes(0, 1) == djs_laxisnum((4, 4, 4), iaxis=1)).all()
        assert (l.swapaxes(0, 2) == djs_laxisnum((4, 4, 4), iaxis=2)).all()
        with raises(ValueError):
            foo = djs_laxisnum((4, 4, 4), iaxis=3)
        #
        # More d
        #
        with raises(ValueError):
            foo = djs_laxisnum((4, 4, 4, 4))

    def test_hogg_iau_name(self):
        assert (hogg_iau_name(354.120375, -0.544777778) ==
                'SDSS J233628.89-003241.2')
        assert (hogg_iau_name(354.120375, -0.544777778, prefix='2MASS') ==
                '2MASS J233628.89-003241.2')
        assert (hogg_iau_name(354.120375, -0.544777778, prefix='') ==
                'J233628.89-003241.2')
        assert (hogg_iau_name(354.120375, -0.544777778, precision=0) ==
                'SDSS J233628.8-003241')
        assert (hogg_iau_name(354.120375, -0.544777778, precision=2) ==
                'SDSS J233628.890-003241.20')
        ra = np.array([354.120375, 7.89439, 36.31915, 110.44730])
        dec = np.array([-0.544777778, -0.35157, 0.47505, 39.35352])
        names = hogg_iau_name(ra, dec)
        assert tuple(names) == ('SDSS J233628.89-003241.2',
                                'SDSS J003134.65-002105.6',
                                'SDSS J022516.59+002830.1',
                                'SDSS J072147.35+392112.6')

    def test_struct_print(self):
        slist = np.zeros((5,), dtype=[('a', 'c16'), ('b', np.bool)])
        with raises(PydlutilsException):
            lines, css = struct_print(slist, silent=True)
        slist = np.array([(1, 2.34, 'five'),
                          (2, 3.456, 'seven'),
                          (3, -4.5678, 'nine')],
                          dtype=[('a', 'i4'), ('bb', 'f4'), ('ccc', 'S5')])
        lines, css = struct_print(slist, silent=True)
        assert lines[0] == 'a bb           ccc  '
        assert lines[1] == '- ------------ -----'
        assert lines[2] == '1         2.34 five '
        assert lines[3] == '2        3.456 seven'
        assert lines[4] == '3      -4.5678 nine '
        assert len(css) == 0
        lines, css = struct_print(slist, silent=True, alias={'ccc': 'c'})
        assert lines[0] == 'a bb           c    '
        assert lines[1] == '- ------------ -----'
        assert lines[2] == '1         2.34 five '
        assert lines[3] == '2        3.456 seven'
        assert lines[4] == '3      -4.5678 nine '
        assert len(css) == 0
        lines, css = struct_print(slist, silent=True,
                                  formatcodes={'a': '{0:02d}'})
        assert lines[0] == 'a  bb           ccc  '
        assert lines[1] == '-- ------------ -----'
        assert lines[2] == '01         2.34 five '
        assert lines[3] == '02        3.456 seven'
        assert lines[4] == '03      -4.5678 nine '
        assert len(css) == 0
        lines, css = struct_print(slist, silent=True, fdigit=3)
        assert lines[0] == 'a bb         ccc  '
        assert lines[1] == '- ---------- -----'
        assert lines[2] == '1       2.34 five '
        assert lines[3] == '2       3.46 seven'
        assert lines[4] == '3      -4.57 nine '
        assert len(css) == 0
        lines, css = struct_print(slist, silent=True, html=True)
        assert lines[0] == '<table>'
        assert lines[1] == '<tr><th>a</th><th>bb</th><th>ccc</th></tr>'
        assert lines[2] == '<tr><td>1</td><td>        2.34</td><td>five </td></tr>'
        assert lines[3] == '<tr><td>2</td><td>       3.456</td><td>seven</td></tr>'
        assert lines[4] == '<tr><td>3</td><td>     -4.5678</td><td>nine </td></tr>'
        assert lines[5] == '</table>'
        assert css[0] == '<style type="text/css">'
        assert css[1] == 'table {'
        assert css[2] == '    border-collapse: collapse;'
        assert css[3] == '}'
        assert css[4] == 'th {'
        assert css[5] == '    padding: 2px;'
        assert css[6] == '    text-align: right;'
        assert css[7] == '    border: 1px solid black;'
        assert css[8] == '    font-weight: bold;'
        assert css[9] == '}'
        assert css[10] == 'td {'
        assert css[11] == '    padding: 2px;'
        assert css[12] == '    text-align: right;'
        assert css[13] == '    border: 1px solid black;'
        assert css[14] == '}'
        assert css[15] == '</style>'
        slist = np.array([(1, 2.34, 'five'),
                          (2, 3.456, 'seven'),
                          (3, -4.5678, 'nine')],
                          dtype=[('a', 'i4'), ('bb', 'f8'), ('ccc', 'S5')])
        lines, css = struct_print(slist, silent=True, ddigit=3)
        assert lines[0] == 'a bb         ccc  '
        assert lines[1] == '- ---------- -----'
        assert lines[2] == '1       2.34 five '
        assert lines[3] == '2       3.46 seven'
        assert lines[4] == '3      -4.57 nine '
        assert len(css) == 0
        with tempfile.NamedTemporaryFile(delete=False) as spf1:
            spf1_name = spf1.name
            lines, css = struct_print(slist, silent=True,
                filename=spf1_name)
        with open(spf1_name, 'rb') as f:
            data = f.read().decode('utf-8')
        assert "\n".join(lines)+"\n" == data
        remove(spf1_name)
        with tempfile.TemporaryFile() as spf2:
            lines, css = struct_print(slist, silent=True, filename=spf2)
            spf2.seek(0)
            data = spf2.read().decode('utf-8')
        assert "\n".join(lines)+"\n" == data
