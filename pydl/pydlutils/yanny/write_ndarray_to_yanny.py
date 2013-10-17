# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def write_ndarray_to_yanny(filename,datatable,structname='mystruct',enums=None,hdr=None):
    """Converts a NumPy record array into a new FTCL/yanny file.

    Returns a new yanny object corresponding to the file.

    Parameters
    ----------
    filename : str
        The name of a parameter file.
    datatable : numpy.ndarray
        A NumPy record array containing data that can be copied into a yanny object.
    structname : str, optional
        The name to give the structure in the yanny file.  Defaults to 'MYSTRUCT'.
    enums : dict, optional
        A dictionary containing enum information.  See details above.
    hdr : dict, optional
        A dictionary containing keyword/value pairs for the 'header' of the yanny file.

    Returns
    -------
    par : yanny
        The yanny object resulting from writing the file.

    Examples
    --------
    """
    from . import yanny
    par = yanny(filename,np=True)
    par['symbols'] = par.dtype_to_struct(datatable.dtype,structname=structname,enums=enums)
    par[structname.upper()] = datatable
    if hdr is not None:
        for key in hdr:
            par[key] = hdr[key]
    par.write(filename)
    return par
