# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy
def write_ndarray_to_yanny(filename,datatables,structnames=None,enums=dict(),hdr=dict()):
    """Converts a NumPy record array into a new FTCL/yanny file.

    Returns a new yanny object corresponding to the file.

    Parameters
    ----------
    filename : str
        The name of a parameter file.
    datatables : numpy.ndarray or tup of numpy.ndarray 
        A NumPy record array containing data that can be copied into a yanny object. A tuple allows
        the writing of multiple data structures to the same file.
    structnames : str or tup of str, optional
        The name to give the structure in the yanny file. If a tuple is given it must be the same
        length as datatables. Defaults to 'MYSTRUCT'.
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

    # Make everything tuples and check that the inputs are correct
    if type(datatables) == numpy.ndarray:
        datatables = (datatables,)
    if not structnames:
        structnames = ['MYSTRUCT{:n}'.format(i) for i in range(len(datatables))]
    if len(datatables) != len(structnames):
        print "The number of data tables is not the same as the number of structnames!"
        return
    
    par = yanny(filename,np=True,debug=True)

    # Now generate the correct info for each structure and add it to the yanny object
    for data, name in zip(datatables, structnames):

        struct = par.dtype_to_struct(data.dtype,structname=name,enums=enums)
        par['symbols']['struct'] += struct['struct']
        par['symbols']['enum'] += struct['enum'] #Yes, this is redundant to do for every struct
        par['symbols'][name.upper()] = struct[name.upper()]
        par[name.upper()] = data

    # write the key/value pairs to the header
    for key in hdr:
        par[key] = hdr[key]
    par.write(filename)

    return par
