# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import os
import numpy as np
from .. import yanny
from .. import write_ndarray_to_yanny
from . import YannyTestCase

class WriteNdarrayTestCase(YannyTestCase):
    """Test class for ndarray to yanny conversions."""
    def test_write_ndarray_to_yanny(self):
        """Test the write_ndarray_to_yanny function."""
        mystruct = np.zeros((4,),
            dtype=[('ra','f8'), ('dec','f8'), ('mag','f4',(5,)),
                ('flags','i4'), ('new_flag','|S5')])
        mystruct['ra'] = np.array([10.0, 20.5, 30.75, 40.55],dtype=np.float64)
        mystruct['dec'] = np.array([-5.1234, -10.74832, 67.994523, 11.437281],dtype=np.float64)
        mystruct['mag'] = np.array([[0.0,1.0,2.0,3.0,4.0],
            [5.1,6.2,7.3,8.4,9.5],
            [22.123,23.95,22.6657,21.0286,22.9876],
            [13.54126,15.37456,14.52647,12.648640,12.0218]],dtype=np.float32)
        mystruct['flags'] = np.array([2**2, 2**4, 2**6, 2**8 + 2**3],dtype=np.int32)
        mystruct['new_flag'] = np.array(['FALSE','TRUE','TRUE','FALSE'],dtype='|S5')
        enums = {'new_flag':('BOOLEAN',('FALSE','TRUE'))}
        par = write_ndarray_to_yanny(self.temp('tempfile1.par'),mystruct,enums=enums)
        assert par['symbols']['enum'][0] == 'typedef enum {\n    FALSE,\n    TRUE\n} BOOLEAN;'
        assert par['symbols']['struct'][0] == 'typedef struct {\n    double ra;\n    double dec;\n    float mag[5];\n    int flags;\n    BOOLEAN new_flag;\n} MYSTRUCT;'
        for k,f in enumerate(('FALSE','TRUE','TRUE','FALSE')):
            assert par['MYSTRUCT']['new_flag'][k] == f
