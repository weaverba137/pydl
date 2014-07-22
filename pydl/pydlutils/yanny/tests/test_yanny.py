# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_yanny():
    """Used to test the yanny class.
    """
    import os
    from .. import yanny
    from numpy import dtype
    #from pydl.pydlutils.yanny import yanny
    #
    # Describe what should be in the object
    #
    pair_dict = {'mjd':'54579','alpha':'beta gamma delta','semicolon':'This pair contains a semicolon;'}
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
    par = yanny(os.path.join(os.path.dirname(__file__),'t','test.par'),np=True)
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
    assert set(par) == set(list(pair_dict) + list(struct_dict) + ['symbols'])
    assert set(par['symbols']) == set(list(struct_dict) + ['struct','enum'])
    assert set(par['symbols']['struct']) == set(symbols)
    assert set(par['symbols']['enum']) == set(enum)
    assert set(par.tables()) == set(struct_dict)
    for t in par.tables():
        assert par.dtype(t) == dtype(struct_dict[t]['dtype'])
        assert par.size(t) == struct_dict[t]['size']
        assert set(par.columns(t)) == set(struct_dict[t]['columns'])
        for c in par.columns(t):
            assert par.type(t,c) == struct_dict[t]['columns'][c]
            #print(par[t][c])
    assert par.isenum('MYSTRUCT','new_flag')
    assert par._enum_cache['BOOLEAN'] == ['FALSE','TRUE']
    assert par._enum_cache['STATUS'] == ['FAILURE','INCOMPLETE','SUCCESS']
    #par.write() # This should fail, since test.par already exists.
    #datatable = {'status_update': {'state':['SUCCESS', 'SUCCESS'],
    #    'timestamp':['2008-06-22 01:27:33','2008-06-22 01:27:36']},
    #    'new_keyword':'new_value'}
    #par.filename = os.path.join(os.path.dirname(__file__),'t','test_append.par')
    #par.append(datatable) # This should also fail, because test_append.par does not exist
    return
#
# Testing purposes
#
if __name__ == '__main__':
    test_yanny()

