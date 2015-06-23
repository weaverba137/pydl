# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_djs_median():
    from numpy import array, median
    from numpy import allclose
    from numpy.random import random, seed
    from astropy.tests.helper import raises
    from .. import djs_median
    seed(424242)
    data = 100.0*random(100)
    data2 = 100.0*random((10,10))
    data3 = 100.0*random((10,10,10))
    data_width_5 = array([ 49.68311576,  74.83671757,  49.68311576,  42.12573137,
        32.11576752,  27.22092569,  15.08905903,  15.08905903,
        27.22092569,  20.07809411,  41.55978953,  41.55978953,
        35.01184343,  35.01184343,  38.7967727 ,  38.7967727 ,
        38.7967727 ,  38.7967727 ,  38.7967727 ,  28.85780425,
        28.85780425,  30.23801035,  30.23801035,  27.0687078 ,
        30.23801035,  64.77058052,  29.42407045,  38.51749618,
        62.48793217,  38.51749618,  29.42407045,  38.51749618,
        18.46919414,  18.46919414,  42.72395431,  54.11584807,
        54.11584807,  75.67593452,  75.67593452,  72.4380356 ,
        61.68761955,  61.68761955,  45.36796557,  45.36796557,
        61.68761955,  76.57621138,  76.57621138,  76.57621138,
        23.28589488,  23.28589488,  13.82755909,  12.10607597,
        13.82755909,  27.51891089,  44.21266068,  44.21266068,
        44.21266068,  47.18083025,  47.18083025,  47.18083025,
        47.18083025,  47.18083025,  35.50809933,  35.50809933,
        25.52222293,  25.52222293,  67.8568479 ,  88.54983822,
        67.8568479 ,  93.40053148,  93.40053148,  64.12514945,
        47.82074715,  47.82074715,  47.82074715,  34.82234113,
        34.82234113,  52.91092248,  78.51075522,  92.16442338,
        92.16442338,  92.16442338,  72.07989558,  72.07989558,
        68.72579501,  72.07989558,  72.07989558,  70.43134908,
        34.55273356,  62.09010468,  62.09010468,  70.43134908,
        68.89705132,  68.89705132,  68.89705132,  66.30426084,
        55.92748086,  55.92748086,  55.92748086,  65.11387444])
    data_width_5_reflect = array([ 49.68311576,  49.68311576,  49.68311576,  42.12573137,
        32.11576752,  27.22092569,  15.08905903,  15.08905903,
        27.22092569,  20.07809411,  41.55978953,  41.55978953,
        35.01184343,  35.01184343,  38.7967727 ,  38.7967727 ,
        38.7967727 ,  38.7967727 ,  38.7967727 ,  28.85780425,
        28.85780425,  30.23801035,  30.23801035,  27.0687078 ,
        30.23801035,  64.77058052,  29.42407045,  38.51749618,
        62.48793217,  38.51749618,  29.42407045,  38.51749618,
        18.46919414,  18.46919414,  42.72395431,  54.11584807,
        54.11584807,  75.67593452,  75.67593452,  72.4380356 ,
        61.68761955,  61.68761955,  45.36796557,  45.36796557,
        61.68761955,  76.57621138,  76.57621138,  76.57621138,
        23.28589488,  23.28589488,  13.82755909,  12.10607597,
        13.82755909,  27.51891089,  44.21266068,  44.21266068,
        44.21266068,  47.18083025,  47.18083025,  47.18083025,
        47.18083025,  47.18083025,  35.50809933,  35.50809933,
        25.52222293,  25.52222293,  67.8568479 ,  88.54983822,
        67.8568479 ,  93.40053148,  93.40053148,  64.12514945,
        47.82074715,  47.82074715,  47.82074715,  34.82234113,
        34.82234113,  52.91092248,  78.51075522,  92.16442338,
        92.16442338,  92.16442338,  72.07989558,  72.07989558,
        68.72579501,  72.07989558,  72.07989558,  70.43134908,
        34.55273356,  62.09010468,  62.09010468,  70.43134908,
        68.89705132,  68.89705132,  68.89705132,  55.92748086,
        65.11387444,  55.92748086,  55.92748086,  55.92748086])
    #
    # Degenerate cases that fall back on numpy.median().
    #
    assert allclose(median(data),djs_median(data))
    assert allclose(median(data2,axis=0),djs_median(data2,dimension=0))
    assert allclose(median(data2,axis=1),djs_median(data2,dimension=1))
    #
    # Test widths.
    #
    assert allclose(data,djs_median(data,width=1))
    # assert_allclose(data_width_5,djs_median(data,width=5))
    # assert_allclose(data_width_5_reflect,djs_median(data,width=5,boundary='reflect'))
    #
    # Exceptions
    #
    with raises(ValueError):
        foo = djs_median(data,width=5,boundary='nearest')
    with raises(ValueError):
        foo = djs_median(data,width=5,boundary='wrap')
    with raises(ValueError):
        foo = djs_median(data,width=5,boundary='foobar')
    with raises(ValueError):
        foo = djs_median(data2,width=5,dimension=1)
    with raises(ValueError):
        foo = djs_median(data3,width=5)
    with raises(ValueError):
        foo = djs_median(data3,width=5,boundary='reflect')
