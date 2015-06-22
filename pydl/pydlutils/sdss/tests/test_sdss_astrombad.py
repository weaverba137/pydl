# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
#
import numpy as np
import pydl.pydlutils.sdss
from astropy.tests.helper import remote_data, raises
#
# First few data lines of opBadfields.par for reference:
#
# BADFIELDS   77   focus  41  60    "Large focus offsets from -74 to 118 microns"
# BADFIELDS   77   focus 106 404    "Large focus offsets from -31 to 165 microns"
# BADFIELDS   77  astrom  30  73     "Large astrometric offset at field 39... 72"
# BADFIELDS   85  astrom   8  28     "Large astrometric offset at field 11... 27"
# BADFIELDS   85 rotator 242 253          "Large rotator offset at field 251 252"
# BADFIELDS  209  astrom   8 116      "Tel. offsets before r-band field 115 -DJS"
# BADFIELDS  209  astrom 137 175       "Tel. offsets after r-band field 145 -DJS"
# BADFIELDS  240 shutter 166 179            "parallel voltages changed frame 176"
# BADFIELDS  250  astrom 456 468               "Large astrometric offset -Manual"
# BADFIELDS  251 shutter  84 103                  "lights on in the dome -Manual"
# BADFIELDS  251   focus  97 160                   "Manual focus tests (donuts!)"
# BADFIELDS  251  astrom 147 159      "Large astrometric offset at field 156 158"
#
class TestSDSSAstrombad(object):
    def setup(self):
        self.opbadfields = np.array([
            (77, 'astrom', 30, 73, 'Large astrometric offset at field 39... 72'),
            (85, 'astrom', 8, 28, 'Large astrometric offset at field 11... 27'),
            (85, 'rotator', 242, 253, 'Large rotator offset at field 251 252'),
            (209, 'astrom', 8, 116, 'Tel. offsets before r-band field 115 -DJS'),
            (209, 'astrom', 137, 175, 'Tel. offsets after r-band field 145 -DJS'),
            (250, 'astrom', 456, 468, 'Large astrometric offset -Manual'),
            (251, 'astrom', 147, 159, 'Large astrometric offset at field 156 158')],
            dtype=[('run', '<i4'), ('problem', 'S8'), ('firstfield', '<i4'), ('lastfield', '<i4'), ('comments', 'S47')])
        pydl.pydlutils.sdss.opbadfields = self.opbadfields

    def teardown(self):
        pass

    def test_sdss_astrombad(self):
        from .. import sdss_astrombad
        assert sdss_astrombad(77,1,20) == False
        assert sdss_astrombad(77,3,35) == True
        assert sdss_astrombad(77,6,77) == False
        assert sdss_astrombad(85,1,15) == True
        assert (sdss_astrombad(np.array([77,85,251]),np.array([1,2,3]),np.array([20,15,151])) == np.array([False,True,True])).all()

    def test_sdss_astrombad_raises(self):
        from .. import sdss_astrombad
        with raises(ValueError):
            foo = sdss_astrombad(77,32,20)
        with raises(ValueError):
            foo = sdss_astrombad(-1,1,20)
        with raises(ValueError):
            foo = sdss_astrombad(2**17,1,20)
        with raises(ValueError):
            foo = sdss_astrombad(-2,1,20)
        with raises(ValueError):
            foo = sdss_astrombad(251,1,2**16)
        with raises(ValueError):
            foo = sdss_astrombad(np.array([77,85,251]),np.array([1]),np.array([20,15,151]))
        with raises(ValueError):
            foo = sdss_astrombad(np.array([77,85,251]),np.array([1,2,3]),np.array([20]))

    @remote_data
    def test_sdss_astrombad_remote(self):
        pydl.pydlutils.sdss.opbadfields = None
        from .. import sdss_astrombad
        assert sdss_astrombad(77,1,20) == False
        assert sdss_astrombad(77,3,35) == True
        assert sdss_astrombad(77,6,77) == False
