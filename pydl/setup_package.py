# Licensed under a 3-clause BSD style license
# -*- coding: utf-8 -*-


def get_package_data():
    # Installs the testing data files.  Unable to get package_data
    # to deal with a directory hierarchy of files, so just explicitly list.
    return {
        'pydl.tests': ['coveragerc', 't/*'],
        'pydl.pydlspec2d.spec2d.tests': ['t/*'],
        'pydl.pydlutils': ['data/filters/*'],
        'pydl.pydlutils.sdss.tests': ['t/*'],
        'pydl.pydlutils.tests': ['t/*'],
        'pydl.pydlutils.trace.tests': ['t/*']
    }


# def requires_2to3():
#     return False
