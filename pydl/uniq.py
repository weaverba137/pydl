#
# $Id$
#
def uniq(x,index=None):
    """Replicates the IDL uniq function.
    """
    import numpy as np
    if index is None:
        indicies = (x != np.roll(x,-1)).nonzero()[0]
        if indicies.size > 0:
            return indicies
        else:
            return np.array([x.size - 1,])
    else:
        q = x[index]
        indicies = (q != np.roll(q,-1)).nonzero()[0]
        if indicies.size > 0:
            return index[indicies]
        else:
            return np.array([q.size - 1,],dtype=index.dtype)

