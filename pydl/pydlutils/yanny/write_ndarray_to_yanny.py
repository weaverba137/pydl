# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def write_ndarray_to_yanny(filename,datatables,structnames=None,
                           enums=None,hdr=None,comments=None):
    """Converts a NumPy record array into a new FTCL/yanny file.

    Returns a new yanny object corresponding to the file.

    Parameters
    ----------
    filename : :class:`str`
        The name of a parameter file.
    datatables : :class:`numpy.ndarray`, :class:`numpy.recarray` or :class:`list` of these.
        A NumPy record array containing data that can be copied into a `yanny` object.
    structnames : :class:`str` or :class:`list` of :class:`str`, optional
        The name(s) to give the structure(s) in the yanny file.  Defaults to 'MYSTRUCT0'.
    enums : :class:`dict`, optional
        A dictionary containing enum information.  See the documentation for
        the :meth:`~pydl.pydlutils.yanny.yanny.dtype_to_struct` method of the yanny object.
    hdr : :class:`dict`, optional
        A dictionary containing keyword/value pairs for the 'header' of the yanny file.
    comments : :class:`str` or :class:`list` of :class:`str`, optional
        A string containing comments that will be added to the start of the new file.

    Returns
    -------
    par : `yanny`
        The `yanny` object resulting from writing the file.

    Raises
    ------
    PydlutilsException
        If `filename` already exists, or if the metadata are incorrect.
    """
    from numpy import ndarray, recarray
    from astropy.extern.six import string_types
    from . import yanny
    from .. import PydlutilsException
    par = yanny(filename,np=True)
    if par:
        #
        # If the file already exists
        #
        raise PydlutilsException("Apparently {0} already exists.".format(filename))
    if isinstance(datatables, (ndarray,recarray)):
        datatables = (datatables,)
    if structnames is None:
        structnames = ["MYSTRUCT{0:d}".format(k) for k in range(len(datatables))]
    if isinstance(structnames, string_types):
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
