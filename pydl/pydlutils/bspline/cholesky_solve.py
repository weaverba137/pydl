#
#
#
def cholesky_solve(a,bb):
    import numpy as np
    b = bb.copy()
    bw = a.shape[0]
    n = b.shape[0] - bw
    kn = bw -1
    spot = np.arange(kn,dtype='i4') + 1
    for j in range(n):
        b[j] /= a[0,j]
        b[j+spot] -= b[j]*a[spot,j]
    spot = kn - np.arange(kn,dtype='i4')
    for j in range(n-1,-1,-1):
        b[j] = (b[j] - np.sum(a[spot,j] * b[j+spot]))/a[0,j]
    return (-1,b)

