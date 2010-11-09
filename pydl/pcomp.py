#
#
#
def pcomp(x,**kwargs):
    """Replacement for IDL pcomp built-in.
    """
    if x.ndim != 2:
        raise ValueError('Input array must be two-dimensional')
    no,nv = x.shape
    if 'standardize' in kwargs:
        xstd = x - np.tile(x.mean(0),no).reshape(x.shape)
        s = np.tile(xstd.std(0),no).reshape(x.shape)
        array = xstd/s
    else:
        array = x
    if 'covariance' in kwargs:
        c = np.cov(array,rowvar=0)
    else:
        c = np.corrcoef(array,rowvar=0)
    #
    # eigh is used for symmetric matrices
    #
    evals, evecs = eigh(c)
    #
    # Sort eigenvalues in descending order
    #
    ie = evals.argsort()[::-1]
    evals = evals[ie]
    evecs = evecs[:,ie]
    #
    # If necessary, add code to fix the signs of the eigenvectors.
    # http://www3.interscience.wiley.com/journal/117912150/abstract
    #
    normevecs = evecs * np.tile(np.sqrt(evals),nv).reshape(nv,nv)
    variances = evals/c.trace()
    derived_data = np.dot(array,normevecs)
    if 'standardize' in kwargs:
        derived_data += xstd
    return {'derived':derived_data,'coefficients':normevecs,
        'variance':variances,'eigenvalues':evals}

