# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import warnings
import numpy as np
from astropy.tests.helper import raises
from astropy.extern import six
from os import chmod
from os.path import dirname, exists, join
from shutil import copy, rmtree
from tempfile import mkdtemp
from time import sleep
from .. import PydlutilsException
from ..yanny import write_ndarray_to_yanny, yanny


class YannyTestCase(object):
    """Based on astropy.io.fits.tests.FitsTestCase.
    """
    save_temp = False

    def setup(self):
        self.data_dir = join(dirname(__file__), 't')
        self.temp_dir = mkdtemp(prefix='yanny-test-')
        # Ignore deprecation warnings--this only affects Python 2.5 and 2.6,
        # since deprecation warnings are ignored by defualt on 2.7
        warnings.simplefilter('ignore')
        warnings.simplefilter('always', UserWarning)
        # raise ValueError("I am setting up a subclass of YannyTestCase!")

    def teardown(self):
        warnings.resetwarnings()
        if not self.save_temp:
            if hasattr(self, 'temp_dir') and exists(self.temp_dir):
                tries = 3
                while tries:
                    try:
                        rmtree(self.temp_dir)
                        break
                    except OSError:
                        # Probably couldn't delete the file because for whatever
                        # reason a handle to it is still open/hasn't been
                        # garbage-collected
                        sleep(0.5)
                        tries -= 1
        # raise ValueError("I am tearing down up a subclass of YannyTestCase!")

    def copy_file(self, filename):
        """Copies a backup of a test data file to the temp dir and sets its
        mode to writeable.
        """
        copy(self.data(filename), self.temp(filename))
        chmod(self.temp(filename), stat.S_IREAD | stat.S_IWRITE)

    def data(self, filename):
        """Returns the path to a test data file."""
        return join(self.data_dir, filename)

    def temp(self, filename):
        """ Returns the full path to a file in the test temp dir."""
        return join(self.temp_dir, filename)


class TestWriteNdarray(YannyTestCase):
    """Test class for ndarray to yanny conversions.
    """
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
        assert par._symbols['enum'][0] == 'typedef enum {\n    FALSE,\n    TRUE\n} BOOLEAN;'
        assert par._symbols['struct'][0] == 'typedef struct {\n    double ra;\n    double dec;\n    float mag[5];\n    int flags;\n    BOOLEAN new_flag;\n} MAGNITUDES;'
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
        assert 'typedef enum {\n    FALSE,\n    TRUE\n} BOOLEAN;' in par._symbols['enum']
        assert 'typedef enum {\n    FAILURE,\n    INCOMPLETE,\n    SUCCESS\n} STATUS;' in par._symbols['enum']
        assert 'typedef struct {\n    double ra;\n    double dec;\n    float mag[5];\n    int flags;\n    BOOLEAN new_flag;\n} MAGNITUDES;' in par._symbols['struct']
        assert 'typedef struct {\n    long timestamp;\n    STATUS state;\n} MY_STATUS;' in par._symbols['struct']
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
            par = write_ndarray_to_yanny(self.data('test.par'), mystruct,
                                        structnames='magnitudes', enums=enums)
        with raises(PydlutilsException):
            par = write_ndarray_to_yanny(self.temp('tempfile3.par'), mystruct,
                                    structnames=('magnitudes', 'my_status'),
                                    enums=enums)

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

    def test_yanny(self):
        """Used to test the yanny class.
        """
        #
        # Describe what should be in the object
        #
        pair_dict = {'mjd': '54579', 'alpha': 'beta gamma delta',
                    'semicolon': 'This pair contains a semicolon;'}
        struct_dict = {
            'MYSTRUCT':{
                'dtype':[('mag', '<f4', (5,)), ('b', 'S33', (5,)), ('foo', 'S25'), ('c', '<f8'), ('flags', '<i4', (2,)), ('new_flag', 'S5')],
                'size':7,
                'columns':{'mag':'float[5]','b':'char[5][]','foo':'char[25]','c':'double','flags':'int[2]','new_flag':'BOOLEAN'},
                },
            'OLD':{
                'dtype':[('foo', '<f4', (3,)), ('bar', 'S10')],
                'size':2,
                'columns':{'foo':'float[3]','bar':'char[10]'},
                },
            'STATUS_UPDATE':{
                'dtype':[('state', 'S10'), ('timestamp', 'S19')],
                'size':11,
                'columns':{'state':'STATUS','timestamp':'char[]'},
                },
            }
        symbols = [
"""typedef struct {
    float mag[5];
    char b[5][];
    char foo[25];
    double c;
    int flags[2];
    BOOLEAN new_flag;
} MYSTRUCT;""",
"""typedef struct {
    float foo<3>; # This is archaic array notation, strongly deprecated,
    char bar<10>; # but still technically supported.
} OLD;""",
"""typedef struct {
    STATUS state;
    char timestamp[]; #UTC timestamp in format 2008-06-21T00:27:33
} STATUS_UPDATE;""",
            ]
        enum = [
"""typedef enum {
    FALSE,
    TRUE
} BOOLEAN;""",
"""typedef enum {
    FAILURE,
    INCOMPLETE,
    SUCCESS
} STATUS;""",
            ]
        #
        # Open the object
        #
        par = yanny(self.data('test.par'))
        #
        # Test the pairs
        #
        assert set(par.pairs()) == set(pair_dict)
        #
        # Test the pair values
        #
        for p in par.pairs():
            assert par[p] == pair_dict[p]
        #
        # Test the structure of the object
        #
        assert set(par) == set(list(pair_dict) + list(struct_dict))
        assert set(par._symbols) == set(list(struct_dict) + ['struct','enum'])
        assert set(par._symbols['struct']) == set(symbols)
        assert set(par._symbols['enum']) == set(enum)
        assert set(par.tables()) == set(struct_dict)
        for t in par.tables():
            assert par.dtype(t) == np.dtype(struct_dict[t]['dtype'])
            assert par.size(t) == struct_dict[t]['size']
            assert set(par.columns(t)) == set(struct_dict[t]['columns'])
            for c in par.columns(t):
                assert par.type(t,c) == struct_dict[t]['columns'][c]
                # print(par[t][c])
        assert par.isenum('MYSTRUCT','new_flag')
        assert par._enum_cache['BOOLEAN'] == ['FALSE','TRUE']
        assert par._enum_cache['STATUS'] == ['FAILURE','INCOMPLETE','SUCCESS']
        #
        # Test values
        #
        assert np.allclose(par['MYSTRUCT'].mag[0],
                            np.array([17.5, 17.546, 17.4, 16.1, 16.0]))
        assert np.allclose(par['MYSTRUCT'].mag[5],
                            np.array([19.3, 18.2, 17.1, 16.0, 15.9]))
        assert par['MYSTRUCT'].foo[1] == six.b("My dog has no nose.")
        assert np.allclose(par['MYSTRUCT'].c[2], 7.24345567)
        assert (par['MYSTRUCT']['flags'][2] == np.array([123123, 0])).all()
        #
        # Test expected write failures.
        #
        # This should fail, since test.par already exists.
        with raises(PydlutilsException):
            par.write()
        datatable = {'status_update': {'state': ['SUCCESS', 'SUCCESS'],
            'timestamp': ['2008-06-22 01:27:33', '2008-06-22 01:27:36']},
            'new_keyword': 'new_value'}
        par.filename = self.data('test_append.par')
        # This should also fail, because test_append.par does not exist.
        with raises(PydlutilsException):
            par.append(datatable)
        return
