# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from astropy.extern import six


def write_ndarray_to_yanny(filename,datatables,structnames=None,
                           enums=None,hdr=None,comments=None):
    """Converts a NumPy record array into a new FTCL/yanny file.

    Returns a new yanny object corresponding to the file.

    Parameters
    ----------
    filename : str
        The name of a parameter file.
    datatables : numpy.ndarray or list of numpy.ndarray
        A NumPy record array containing data that can be copied into a yanny object.
    structnames : str or list of str, optional
        The name(s) to give the structure(s) in the yanny file.  Defaults to 'MYSTRUCT0'.
    enums : dict, optional
        A dictionary containing enum information.  See the documentation for
        the `dtype_to_struct` method of the yanny object.
    hdr : dict, optional
        A dictionary containing keyword/value pairs for the 'header' of the yanny file.
    comments : str or list of str, optional
        A string containing comments that will be added to the start of the new file.

    Returns
    -------
    par : pydl.pydlutils.yanny.yanny
        The `yanny` object resulting from writing the file.

    Examples
    --------
    """
    from numpy import ndarray
    from . import yanny
    from .. import PydlutilsException
    par = yanny(filename,np=True)
    if par:
        #
        # If the file already exists
        #
        raise PydlutilsException("Apparently {0} already exists.".format(filename))
    if type(datatables) == ndarray:
        datatables = (datatables,)
    if structnames is None:
        structnames = ["MYSTRUCT{0:d}".format(k) for k in range(len(datatables))]
    if isinstance(structnames, six.string_types):
        structnames = (structnames,)
    if len(datatables) != len(structnames):
        raise PydlutilsException("The data tables and their names do not match!")
    for k in range(len(datatables)):
        struct = par.dtype_to_struct(datatables[k].dtype,structname=structnames[k],enums=enums)
        par['symbols']['struct'] += struct['struct']
        par['symbols'][structnames[k].upper()] = struct[structnames[k].upper()]
        if enums is not None and len(par['symbols']['enum']) == 0:
            par['symbols']['enum'] = struct['enum']
        par[structnames[k].upper()] = datatables[k]
    if hdr is not None:
        for key in hdr:
            par[key] = hdr[key]
    par.write(filename,comments=comments)
    return par
