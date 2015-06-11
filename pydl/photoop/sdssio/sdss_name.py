# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def sdss_name(ftype, run, camcol, field, rerun='', thisfilter='r',no_path=False):
    """Return the name of an SDSS data file including path.

    Parameters
    ----------
    ftype : str
        The general type of the file, for example ``'reObj'``
    run : int
        The run number.
    camcol : int
        The camcol number.
    field : int
        The field number
    rerun : str, optional
        If necessary, set the rerun name using this argument.
    thisfilter : int or str, optional
        If necessary, set the filter using this argument.
    no_path : bool, optional
        Normally, sdss_name returns the full path.  If `no_path` is ``True``, only
        the basename of the file is returned.

    Returns
    -------
    sdss_name : str
        The full file name, normally including the full path.

    Raises
    ------
    KeyError
        If the file type is unknown.

    """
    from os import getenv
    from os.path import join
    from . import _name_formats, sdss_path, filtername
    # cname = ('u','g','r','i','z')
    # camrow = (3,5,1,2,4)
    if ftype == 'reObj':
        if getenv('PHOTO_RESOLVE') is None:
            myftype = 'reObjRun'
        else:
            myftype = 'reObjGlobal'
    else:
        myftype = ftype
    thisfilter = filtername(thisfilter)
    indict = {
        'ftype':myftype, 'run':run, 'camcol':camcol, 'field':field,
        'filter':thisfilter, 'rerun':rerun
        }
    try:
        fullname = _name_formats[myftype].format(**indict)
    except KeyError:
        raise KeyError("Unknown FTYPE = {0}".format(myftype))
    if not no_path:
        datadir = sdss_path(myftype, run, camcol, rerun)
        fullname = join(datadir, fullname)
    return fullname
