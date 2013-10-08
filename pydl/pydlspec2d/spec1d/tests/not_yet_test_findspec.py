#
# $Id$
#
def test_findspec():
    from pydlspec2d.spec1d import findspec
    slist = findspec(infile='file.in',sdss=True)
    if slist is not None:
        print slist
    return
