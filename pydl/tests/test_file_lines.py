# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals
#
def test_file_lines():
    from ..file_lines import file_lines
    from os.path import basename, dirname, join
    import glob
    #
    # Find the test files
    #
    datadir = join(dirname(__file__),'t')
    fileglob = join(datadir,'this-file-contains-*-lines.txt')
    plainfiles = glob.glob(fileglob)
    gzfiles = glob.glob(fileglob+'.gz')
    for p in plainfiles:
        n = file_lines(p)
        number_of_lines = int(basename(p).split('-')[3])
        assert n == number_of_lines
    for p in gzfiles:
        n = file_lines(p,compress=True)
        number_of_lines = int(basename(p).split('-')[3])
        assert n == number_of_lines
    #
    # Test list passing
    #
    n = file_lines(plainfiles)
    number_of_lines = [int(basename(p).split('-')[3]) for p in plainfiles]
    assert n == number_of_lines
    n = file_lines(gzfiles,compress=True)
    number_of_lines = [int(basename(p).split('-')[3]) for p in gzfiles]
    assert n == number_of_lines
    #
    # Make sure empty files work
    #
    n = file_lines(join(datadir,'this-file-is-empty.txt'))
    assert n == 0
