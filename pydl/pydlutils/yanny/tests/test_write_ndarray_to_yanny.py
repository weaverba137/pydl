# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import os
import numpy as np
from .. import yanny
from .. import write_ndarray_to_yanny
from ... import PydlutilsException
from . import YannyTestCase
from astropy.tests.helper import raises

class TestWriteNdarray(YannyTestCase):
    """Test class for ndarray to yanny conversions."""
    save_temp = False
    def test_write_single_ndarray_to_yanny(self):
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
        par = write_ndarray_to_yanny(self.temp('tempfile1.par'),mystruct,structnames='magnitudes',enums=enums)
        assert par['symbols']['enum'][0] == 'typedef enum {\n    FALSE,\n    TRUE\n} BOOLEAN;'
        assert par['symbols']['struct'][0] == 'typedef struct {\n    double ra;\n    double dec;\n    float mag[5];\n    int flags;\n    BOOLEAN new_flag;\n} MAGNITUDES;'
        for k,f in enumerate(('FALSE','TRUE','TRUE','FALSE')):
            assert par['MAGNITUDES']['new_flag'][k].decode() == f

    def test_write_multiple_ndarray_to_yanny(self):
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
        status = np.zeros((4,),dtype=[('timestamp','i8'),('state','S10')])
        status['timestamp'] = np.array([1382384327,1382384527,1382384727,1382384927],dtype='i8')
        status['state'] = np.array(['SUCCESS','SUCCESS','FAILURE','INCOMPLETE'],dtype='S10')
        enums = {'new_flag':('BOOLEAN',('FALSE','TRUE')),'state':('STATUS',('FAILURE','INCOMPLETE','SUCCESS'))}
        par = write_ndarray_to_yanny(self.temp('tempfile2.par'),(mystruct,status),structnames=('magnitudes','my_status'),enums=enums,comments=['This is a test','This is another test'])
        assert 'typedef enum {\n    FALSE,\n    TRUE\n} BOOLEAN;' in par['symbols']['enum']
        assert 'typedef enum {\n    FAILURE,\n    INCOMPLETE,\n    SUCCESS\n} STATUS;' in par['symbols']['enum']
        assert 'typedef struct {\n    double ra;\n    double dec;\n    float mag[5];\n    int flags;\n    BOOLEAN new_flag;\n} MAGNITUDES;' in par['symbols']['struct']
        assert 'typedef struct {\n    long timestamp;\n    STATUS state;\n} MY_STATUS;' in par['symbols']['struct']
        for k,f in enumerate(('FALSE','TRUE','TRUE','FALSE')):
            assert par['MAGNITUDES']['new_flag'][k].decode() == f

    def test_write_ndarray_to_yanny_exceptions(self):
        """Make sure certain execptions are raised by the write_ndarray_to_yanny function."""
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
        with raises(PydlutilsException):
            par = write_ndarray_to_yanny(os.path.join(self.data_dir,'test.par'),mystruct,structnames='magnitudes',enums=enums)
        with raises(PydlutilsException):
            par = write_ndarray_to_yanny(self.temp('tempfile3.par'),mystruct,structnames=('magnitudes','my_status'),enums=enums)

    def test_write_ndarray_to_yanny_misc(self):
        """Test other minor features of write_ndarray_to_yanny."""
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
        hdr = {'keyword1':'value1','keyword2':'value2'}
        par = write_ndarray_to_yanny(self.temp('tempfile1.par'),mystruct,hdr=hdr,enums=enums)
        assert 'MYSTRUCT0' in par
        assert par['keyword1'] == hdr['keyword1']
        assert par['keyword2'] == hdr['keyword2']
