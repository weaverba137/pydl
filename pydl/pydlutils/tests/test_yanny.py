# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import warnings
import json
from os import chmod, remove
from os.path import exists, join
from shutil import copy, rmtree
from tempfile import mkdtemp
from time import sleep
from collections import OrderedDict
import numpy as np
from astropy.tests.helper import catch_warnings, raises
from astropy.utils.data import get_pkg_data_filename
from astropy.table import Table
from astropy.io.registry import (register_identifier, register_reader,
                                 register_writer)
from .. import PydlutilsException, PydlutilsUserWarning
from ..yanny import (write_ndarray_to_yanny, yanny, is_yanny,
                     read_table_yanny, write_table_yanny)


register_identifier('yanny', Table, is_yanny)
register_reader('yanny', Table, read_table_yanny)
register_writer('yanny', Table, write_table_yanny)


class YannyTestCase(object):
    """Based on astropy.io.fits.tests.FitsTestCase.
    """
    save_temp = False

    def setup(self):
        self.temp_dir = mkdtemp(prefix='yanny-test-')
        # Ignore deprecation warnings--this only affects Python 2.5 and 2.6,
        # since deprecation warnings are ignored by defualt on 2.7
        warnings.simplefilter('ignore')
        warnings.simplefilter('always', UserWarning)
        # raise ValueError("I am setting up a subclass of YannyTestCase!")
        with open(get_pkg_data_filename("t/yanny_data.json")) as js:
            self.test_data = json.load(js)

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
                        # Probably couldn't delete the file because for
                        # whatever reason a handle to it is still open/hasn't
                        # been garbage-collected.
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
        """Returns the path to a test data file.
        """
        return get_pkg_data_filename('t/'+filename)

    def temp(self, filename):
        """Returns the full path to a file in the test temp dir.
        """
        return join(self.temp_dir, filename)

    def json2dtype(self, data):
        """Convert JSON-encoded dtype data into a real dtype object.
        """
        stuff = list()
        for k in data:
            if len(k) == 3:
                stuff.append(tuple([str(k[0]), str(k[1]), tuple(k[2])]))
            else:
                stuff.append(tuple([str(k[0]), str(k[1])]))
        return np.dtype(stuff)


class TestYanny(YannyTestCase):
    """Test class for yanny files.
    """
    save_temp = False

    def test_is_yanny(self):
        """Test the yanny identifier.
        """
        par = yanny(self.data('test.par'))
        assert is_yanny('read', None, None, par)
        assert is_yanny('read', 'test.par', None)
        with open(self.data('test.par'), 'rb') as fileobj:
            assert is_yanny('read', None, fileobj)

    def test_read_table_yanny(self):
        """Test reading to an astropy Table.
        """
        filename = self.data('test_table.par')
        with raises(PydlutilsException):
            t = Table.read(filename)
        with raises(KeyError):
            t = Table.read(filename, tablename='foo')
        t = Table.read(filename, tablename='test')
        assert isinstance(t.meta, OrderedDict)
        assert t.meta['name'] == 'first table'
        assert (t['a'] == np.array([1, 4, 5])).all()
        assert (t['c'] == np.array([b'x', b'y', b'z'])).all()

    def test_write_table_yanny(self):
        """Test writing an astropy Table to yanny.
        """
        filename = self.temp('table.par')
        a = [1, 4, 5]
        b = [2.0, 5.0, 8.2]
        c = [b'x', b'y', b'z']
        t = Table([a, b, c], names=('a', 'b', 'c'),
                  meta={'name': 'first table'},
                  dtype=('i8', 'f8', 'S1'))
        t.write(filename, tablename='test')
        par1 = yanny(filename)
        par2 = yanny(self.data('test_table.par'))
        assert par1 == par2
        remove(filename)

    def test_write_single_ndarray_to_yanny(self):
        """Test the write_ndarray_to_yanny function.
        """
        test_data = self.test_data['tempfile1.par']
        table = test_data['structures']['MYSTRUCT0']
        dt = self.json2dtype(table['dtype'])
        mystruct = np.zeros((table['size'],), dtype=dt)
        for col in table['dtype']:
            mystruct[col[0]] = np.array(table['data'][col[0]], dtype=col[1])
        enums = {'new_flag': test_data['enums']['new_flag']}
        par = write_ndarray_to_yanny(self.temp('tempfile1.par'),
                                     mystruct,
                                     structnames='magnitudes',
                                     enums=enums)
        assert (par._symbols['enum'][0] ==
                'typedef enum {\n    FALSE,\n    TRUE\n} BOOLEAN;')
        assert par._symbols['struct'][0] == 'typedef struct {\n    double ra;\n    double dec;\n    float mag[5];\n    int flags;\n    BOOLEAN new_flag;\n} MAGNITUDES;'
        for k, f in enumerate(('FALSE', 'TRUE', 'TRUE', 'FALSE')):
            assert par['MAGNITUDES']['new_flag'][k].decode() == f

    def test_write_multiple_ndarray_to_yanny(self):
        """Test the write_ndarray_to_yanny function.
        """
        test_data = self.test_data['tempfile1.par']
        table = test_data['structures']['MYSTRUCT0']
        dt = self.json2dtype(table['dtype'])
        mystruct = np.zeros((table['size'],), dtype=dt)
        for col in table['dtype']:
            mystruct[col[0]] = np.array(table['data'][col[0]], dtype=col[1])
        enums = test_data['enums']
        table = test_data['structures']['MY_STATUS']
        dt = self.json2dtype(table['dtype'])
        status = np.zeros((table['size'],), dtype=dt)
        for col in table['dtype']:
            status[col[0]] = np.array(table['data'][col[0]], dtype=col[1])
        enums = test_data['enums']
        par = write_ndarray_to_yanny(self.temp('tempfile2.par'),
                                    (mystruct, status),
                                    structnames=('magnitudes', 'my_status'),
                                    enums=enums,
                                    comments=['This is a test',
                                    'This is another test'])
        assert ('typedef enum {\n    FALSE,\n    TRUE\n} BOOLEAN;' in
                par._symbols['enum'])
        assert 'typedef enum {\n    FAILURE,\n    INCOMPLETE,\n    SUCCESS\n} STATUS;' in par._symbols['enum']
        assert 'typedef struct {\n    double ra;\n    double dec;\n    float mag[5];\n    int flags;\n    BOOLEAN new_flag;\n} MAGNITUDES;' in par._symbols['struct']
        assert 'typedef struct {\n    long timestamp;\n    STATUS state;\n} MY_STATUS;' in par._symbols['struct']
        for k, f in enumerate(('FALSE', 'TRUE', 'TRUE', 'FALSE')):
            assert par['MAGNITUDES']['new_flag'][k].decode() == f
        for k, f in enumerate(["SUCCESS", "SUCCESS", "FAILURE", "INCOMPLETE"]):
            assert par['MY_STATUS']['state'][k].decode() == f

    def test_write_ndarray_to_yanny_exceptions(self):
        """Make sure certain execptions are raised by the
        write_ndarray_to_yanny function.
        """
        test_data = self.test_data['tempfile1.par']
        table = test_data['structures']['MYSTRUCT0']
        dt = self.json2dtype(table['dtype'])
        mystruct = np.zeros((table['size'],), dtype=dt)
        for col in table['dtype']:
            mystruct[col[0]] = np.array(table['data'][col[0]], dtype=col[1])
        enums = {'new_flag': test_data['enums']['new_flag']}
        with raises(PydlutilsException):
            par = write_ndarray_to_yanny(self.data('test.par'), mystruct,
                                        structnames='magnitudes', enums=enums)
        with raises(PydlutilsException):
            par = write_ndarray_to_yanny(self.temp('tempfile3.par'), mystruct,
                                    structnames=('magnitudes', 'my_status'),
                                    enums=enums)

    def test_write_ndarray_to_yanny_misc(self):
        """Test other minor features of write_ndarray_to_yanny.
        """
        test_data = self.test_data["tempfile1.par"]
        table = test_data['structures']['MYSTRUCT0']
        dt = self.json2dtype(table['dtype'])
        mystruct = np.zeros((table['size'],), dtype=dt)
        for col in table['dtype']:
            mystruct[col[0]] = np.array(table['data'][col[0]], dtype=col[1])
        enums = {'new_flag': test_data['enums']['new_flag']}
        hdr = test_data['pairs']
        par = write_ndarray_to_yanny(self.temp('tempfile1.par'), mystruct,
                                    hdr=hdr, enums=enums)
        assert 'MYSTRUCT0' in par
        assert par['keyword1'] == hdr['keyword1']
        assert par['keyword2'] == hdr['keyword2']

    def test_yanny(self):
        """Used to test the yanny class.
        """
        #
        # Describe what should be in the object
        #
        pair_data = self.test_data["test.par"]["pairs"]
        struct_data = self.test_data["test.par"]["structures"]
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
        with open(self.data('test.par')) as f:
            file_data = f.read()
        #
        # Open the object
        #
        par = yanny(self.data('test.par'))
        #
        # Test some static methods, etc.
        #
        assert par.get_token('abcd') == ('abcd', '')
        assert str(par) == file_data
        assert repr(par) == file_data
        par2 = yanny(self.data('test_table.par'))
        assert par != par2
        #
        # Test types
        #
        assert par.type('MYSTRUCT', 'c') == 'double'
        assert par.type('FOOBAR', 'd') is None
        assert par.type('MYSTRUCT', 'foobar') is None
        assert not par2.isenum('TEST', 'a')
        assert par.array_length('MYSTRUCT', 'c') == 1
        assert par.char_length('MYSTRUCT', 'c') is None
        #
        # Test the pairs
        #
        assert set(par.pairs()) == set(pair_data)
        #
        # Test the pair values
        #
        for p in par.pairs():
            assert par[p] == pair_data[p]
        #
        # Test the structure of the object
        #
        assert (set(par) == set(list(pair_data) +
                list(struct_data)))
        assert (set(par._symbols) == set(list(struct_data) +
                ['struct', 'enum']))
        assert set(par._symbols['struct']) == set(symbols)
        assert set(par._symbols['enum']) == set(enum)
        assert set(par.tables()) == set(struct_data)
        for t in par.tables():
            assert (par.dtype(t) ==
                    self.json2dtype(struct_data[t]['dtype']))
            assert par.size(t) == struct_data[t]['size']
            assert (set(par.columns(t)) ==
                    set(struct_data[t]['columns']))
            for c in par.columns(t):
                assert (par.type(t, c) ==
                        struct_data[t]['columns'][c])
                # print(par[t][c])
        assert par.isenum('MYSTRUCT', 'new_flag')
        assert par._enum_cache['BOOLEAN'] == ['FALSE', 'TRUE']
        assert (par._enum_cache['STATUS'] ==
                ['FAILURE', 'INCOMPLETE', 'SUCCESS'])
        #
        # Test values
        #
        assert np.allclose(par['MYSTRUCT'].mag[0],
                            np.array([17.5, 17.546, 17.4, 16.1, 16.0]))
        assert np.allclose(par['MYSTRUCT'].mag[5],
                            np.array([19.3, 18.2, 17.1, 16.0, 15.9]))
        assert par['MYSTRUCT'].foo[1] == b"My dog has no nose."
        assert np.allclose(par['MYSTRUCT'].c[2], 7.24345567)
        assert (par['MYSTRUCT']['flags'][2] == np.array([123123, 0])).all()
        #
        # Test expected write failures.
        #
        # This should fail, since test.par already exists.
        with raises(PydlutilsException):
            par.write()
        with catch_warnings(PydlutilsUserWarning) as w:
            par.append({})
        datatable = {'status_update': {'state': ['SUCCESS', 'SUCCESS'],
            'timestamp': ['2008-06-22 01:27:33', '2008-06-22 01:27:36']},
            'new_keyword': 'new_value'}
        par.filename = self.temp('test_append.par')
        # This should also fail, because test_append.par does not exist.
        with raises(PydlutilsException):
            par.append(datatable)
        return
