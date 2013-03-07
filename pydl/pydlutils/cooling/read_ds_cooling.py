#
#
#
def read_ds_cooling(fname,logT=None):
    """Read in Dopita & Sutherland 1993 cooling function
    """
    import os
    import os.path
    import numpy as np
    from pydlutils import PydlutilsException
    from pydlutils.goddard.misc import readcol
    path = os.path.join(os.getenv('IDLUTILS_DIR'),'data','cooling',fname)
    if os.path.exists(path):
        logTin, nelec, nH, nt, loglambdanet, loglambdanorm = readcol(fname,skip=4)
        if logT is None:
            return (logTin, loglambdanet)
        else:
            loglambda = np.interp(logT,logTin,loglambdanet)
            return (logT,loglambda)
    else:
        raise PydlutilsException("Could not find "+fname)
