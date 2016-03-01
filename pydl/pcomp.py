# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
from astropy.utils import lazyproperty


class pcomp(object):
    """Replicates the IDL ``PCOMP()`` function.

    The attributes of this class are all read-only properties, implemented
    with :class:`~astropy.utils.decorators.lazyproperty`.

    Parameters
    ----------
    x : array-like
        A 2-D array with :math:`N` rows and :math:`M` columns.
    standardize : :class:`bool`, optional
        If set to ``True``, the input data will have its mean subtracted off
        and will be scaled to unit variance.
    covariance : :class:`bool`, optional.
        If set to ``True``, the covariance matrix of the data will be used for
        the computation.  Otherwise the correlation matrix will be used.

    Notes
    -----

    References
    ----------
    http://www.exelisvis.com/docs/PCOMP.html

    Examples
    --------
    """

    def __init__(self, x, standardize=False, covariance=False):
        from scipy.linalg import eigh
        if x.ndim != 2:
            raise ValueError('Input array must be two-dimensional')
        no, nv = x.shape
        self._nv = nv
        if standardize:
            xstd = x - np.tile(x.mean(0), no).reshape(x.shape)
            s = np.tile(xstd.std(0), no).reshape(x.shape)
            self._array = xstd/s
            self._xstd = xstd
        else:
            self._array = x
            self._xstd = None
        self._standardize = standardize
        if covariance:
            self._c = np.cov(self._array, rowvar=0)
        else:
            self._c = np.corrcoef(self._array, rowvar=0)
        self._covariance = covariance
        #
        # eigh is used for symmetric matrices
        #
        evals, evecs = eigh(self._c)
        #
        # Sort eigenvalues in descending order
        #
        ie = evals.argsort()[::-1]
        self._evals = evals[ie]
        self._evecs = evecs[:, ie]
        #
        # If necessary, add code to fix the signs of the eigenvectors.
        # http://www3.interscience.wiley.com/journal/117912150/abstract
        #
        return

    @lazyproperty
    def coefficients(self):
        """(:class:`~numpy.ndarray`) The principal components.
        These are the coefficients of `derived`.
        Basically, they are a re-scaling of the eigenvectors.
        """
        return self._evecs * np.tile(np.sqrt(self._evals), self._nv).reshape(
            self._nv, self._nv)

    @lazyproperty
    def derived(self):
        """(:class:`~numpy.ndarray`) The derived variables.
        """
        derived_data = np.dot(self._array, self.coefficients)
        if self._standardize:
            derived_data += self._xstd
        return derived_data

    @lazyproperty
    def variance(self):
        """(:class:`~numpy.ndarray`) The variances of each derived variable.
        """
        return self._evals/self._c.trace()

    @lazyproperty
    def eigenvalues(self):
        """(:class:`~numpy.ndarray`) The eigenvalues.
        There is one eigenvalue for each principal component.
        """
        return self._evals
