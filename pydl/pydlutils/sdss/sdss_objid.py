# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def sdss_objid(run,camcol,field,objnum,rerun=301,skyversion=None):
    """Convert SDSS photometric identifiers into CAS-style ObjID.

    Bits are assigned in objid thus:

    ===== ========== ===============================================
    Bits  Name       Comment
    ===== ========== ===============================================
    63    empty      unassigned
    59-62 skyVersion resolved sky version (0-15)
    48-58 rerun      number of pipeline rerun
    32-47 run        run number
    29-31 camcol     camera column (1-6)
    28    firstField [is this the first field in segment?] 0 for now
    16-27 field      field number within run
    0-15  object     object number within field
    ===== ========== ===============================================

    Parameters
    ----------
    run, camcol, field, objnum : :class:`int` or array of int
        Run, camcol, field and object number within field.  If arrays are passed,
        all must have the same length.
    rerun, skyversion : :class:`int` or array of int, optional
        Rerun and skyversion usually don't change very much.  If supplied, make
        sure the size matches all the other values.

    Returns
    -------
    sdss_objid : :class:`numpy.ndarray` of :class:`numpy.int64`
        The ObjIDs of the objects.

    Raises
    ------
    ValueError
        If the sizes of the arrays don't match or if the array values are out of bounds.

    Notes
    -----
    firstField flag never set.

    Examples
    --------
    >>> from pydl.pydlutils.sdss import sdss_objid
    >>> sdss_objid(3704,3,91,146)
    array([1237661382772195474])
    """
    from . import default_skyversion
    from numpy import array, int64, zeros
    if skyversion is None:
        skyversion = default_skyversion()
    if isinstance(run,int):
        run = array([run],dtype=int64)
    if isinstance(camcol,int):
        camcol = array([camcol],dtype=int64)
    if isinstance(field,int):
        field = array([field],dtype=int64)
    if isinstance(objnum,int):
        objnum = array([objnum],dtype=int64)
    if isinstance(rerun,int):
        if rerun == 301:
            rerun = zeros(run.shape,dtype=int64) + 301
        else:
            rerun = array([rerun],dtype=int64)
    if isinstance(skyversion,int):
        if skyversion == default_skyversion():
            skyversion = zeros(run.shape,dtype=int64) + default_skyversion()
        else:
            skyversion = array([skyversion],dtype=int64)

    #
    # Check that all inputs have the same shape.
    #
    firstfield = zeros(run.shape,dtype=int64)
    if run.shape != camcol.shape:
        raise ValueError("camcol.shape does not match run.shape!")
    if run.shape != field.shape:
        raise ValueError("field.shape does not match run.shape!")
    if run.shape != objnum.shape:
        raise ValueError("objnum.shape does not match run.shape!")
    if run.shape != rerun.shape:
        raise ValueError("rerun.shape does not match run.shape!")
    if run.shape != skyversion.shape:
        raise ValueError("skyversion.shape does not match run.shape!")
    #
    # Check ranges of parameters
    #
    if ((skyversion < 0) | (skyversion >= 16)).any():
        raise ValueError("skyversion values are out-of-bounds!")
    if ((rerun < 0) | (rerun >= 2**11)).any():
        raise ValueError("rerun values are out-of-bounds!")
    if ((run < 0) | (run >= 2**16)).any():
        raise ValueError("run values are out-of-bounds!")
    if ((camcol < 1) | (camcol > 6)).any():
        raise ValueError("camcol values are out-of-bounds!")
    if ((field < 0) | (field >= 2**12)).any():
        raise ValueError("camcol values are out-of-bounds!")
    if ((objnum < 0) | (objnum >= 2**16)).any():
        raise ValueError("id values are out-of-bounds!")
    #
    # Compute the objid
    #
    objid = ((skyversion << 59) |
        (rerun << 48) |
        (run << 32) |
        (camcol << 29) |
        (firstfield << 28) |
        (field << 16) |
        (objnum))
    return objid
