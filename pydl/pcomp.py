# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
from astropy.utils import lazyproperty
#
class pcomp(object):
    """Replicates the IDL PCOMP() function.

    Parameters
    ----------
    x : array_like
        A 2-D arrray.
    standardize : bool, optional
        If set to ``True``, The input data will have its mean subtracted off
        and will be scaled to uint standard deviation.
    covariance : bool, optional.
        If set to ``True``, the covariance matrix of the data will be used for
        the computation.  Otherwise the correlation matrix will be used.

    Attributes
    ----------
    coefficients
    derived
    variance
    eigenvalues

    Notes
    -----

    References
    ----------
    http://www.exelisvis.com/docs/PCOMP.html

    Examples
    --------
    """
    def __init__(self,x,standardize=False,covariance=False):
        from scipy.linalg import eigh
        if x.ndim != 2:
            raise ValueError('Input array must be two-dimensional')
        no,nv = x.shape
        self._nv = nv
        if standardize:
            xstd = x - np.tile(x.mean(0),no).reshape(x.shape)
            s = np.tile(xstd.std(0),no).reshape(x.shape)
            self._array = xstd/s
            self._xstd = xstd
        else:
            self._array = x
            self._xstd = None
        self._standardize = standardize
        if covariance:
            self._c = np.cov(array,rowvar=0)
        else:
            self._c = np.corrcoef(array,rowvar=0)
        self._covariance = covariance
        #
        # eigh is used for symmetric matrices
        #
        evals, evecs = eigh(c)
        #
        # Sort eigenvalues in descending order
        #
        ie = evals.argsort()[::-1]
        self._evals = evals[ie]
        self._evecs = evecs[:,ie]
        #
        # If necessary, add code to fix the signs of the eigenvectors.
        # http://www3.interscience.wiley.com/journal/117912150/abstract
        #
        return
    #
    #
    #
    @lazyproperty
    def coefficients(self):
        return self._evecs * np.tile(np.sqrt(self._evals),self._nv).reshape(self._nv,self._nv)
    #
    #
    #
    @lazyproperty
    def derived(self):
        derived_data = np.dot(self._array,self.coefficients)
        if self._standardize:
            derived_data += self._xstd
        return derived_data
    #
    #
    #
    @lazyproperty
    def variance(self):
        return self._evals/self._c.trace()
    #
    #
    #
    @lazyproperty
    def eigenvalues(self):
        return self._evals
