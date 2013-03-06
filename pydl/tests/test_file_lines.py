# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_file_lines():
    from ..file_lines import file_lines
    import glob
    #
    # Find the test files
    #
    plainfiles = glob.glob('t/this-file-contains-*-lines.txt')
    gzfiles = glob.glob('t/this-file-contains-*-lines.txt.gz')
    for p in plainfiles:
        n = file_lines(p)
        number_of_lines = int(p.split('-')[3])
        assert n == number_of_lines
    for p in gzfiles:
        n = file_lines(p,compressed=True)
        number_of_lines = int(p.split('-')[3])
        assert n == number_of_lines
    #
    # Test list passing
    #
    n = file_lines(plainfiles)
    number_of_lines = [int(p.split('-')[3]) for p in plainfiles]
    assert n == number_of_lines
    n = file_lines(gzfiles,compress=True)
    number_of_lines = [int(p.split('-')[3]) for p in gzfiles]
    assert n == number_of_lines
