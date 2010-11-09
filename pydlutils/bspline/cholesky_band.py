def cholesky_band(l,mininf=0.0,verbose=False):
    import numpy as np
    lower = l.copy()
    bw,nn = lower.shape
    n = nn - bw
    # print lower[0,0:n]
    negative = lower[0,0:n] <= mininf
    if negative.any() or not np.all(np.isfinite(lower)):
        if verbose:
            print 'Bad entries.'
            print negative.nonzero()[0]
        return (negative.nonzero()[0],l)
    kn = bw - 1
    spot = np.arange(kn,dtype='i4') + 1
    bi = np.arange(kn,dtype='i4')
    for i in range(1,kn):
        bi = np.append(bi, np.arange(kn-i,dtype='i4') + (kn+1)*i)
    for j in range(n):
        lower[0,j] = np.sqrt(lower[0,j])
        lower[spot,j] /= lower[0,j]
        x = lower[spot,j]
        if not np.all(np.isfinite(x)):
            if verbose:
                print 'NaN found in cholesky_band.'
            return (j,l)
        hmm = np.outer(x,x)
        here = bi+(j+1)*bw
        lower.T.flat[here] -= hmm.flat[bi]
    return (-1,lower)

