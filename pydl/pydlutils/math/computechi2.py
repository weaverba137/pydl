# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
#
class computechi2(object):
    """Solve the linear set of equations Ax=b using SVD.

    The attributes of this class are all implemented as read-only properties.

    Parameters
    ----------
    bvec : ndarray
        The b vector in Ax=b. This vector has length N.
    sqivar : ndarray
        The reciprocal of the errors in `b`.  The name comes from the square
        root of the inverse variance, which is what this is.
    amatrix : ndarray
        The matrix A in Ax=b.  The shape of this matrix is (N,M).

    Attributes
    ----------
    chi2 : float
        The chi**2 value of the fit.
    acoeff : ndarray
        The fit parameters, x, in Ax=b.  This vector has length M.
    yfit : ndarray
        The evaluated best-fit at each point.  This vector has length N.
    dof : int
        The degrees of freedom of the fit.  This is the number of values of `bvec`
        that have `sqivar` > 0 minus the number of fit paramaters, which is
        equal to M.
    covar : ndarray
        The covariance matrix.  The shape of this matrix is (M,M).
    var : ndarray
        The variances of the fit.  This is identical to the diagonal of the
        covariance matrix.  This vector has length M.
    """
    #
    #
    #
    def __init__(self,bvec,sqivar,amatrix):
        """Initialize the object and perform initial computations.
        """
        from numpy.linalg import svd
        #
        # Save the inputs
        #
        self.bvec = bvec
        self.sqivar = sqivar
        self.amatrix = amatrix
        if len(amatrix.shape) > 1:
            self.nstar = amatrix.shape[1]
        else:
            self.nstar = 1
        self.bvec = bvec*sqivar
        self.mmatrix = self.amatrix * np.tile(sqivar,self.nstar).reshape(self.nstar,bvec.size).transpose()
        mm = np.dot(self.mmatrix.T,self.mmatrix)
        self.uu,self.ww,self.vv = svd(mm,full_matrices=False)
        self.mmi = np.dot((self.vv.T / np.tile(self.ww,self.nstar).reshape(self.nstar,self.nstar)),self.uu.T)
        return
    #
    #
    #
    @property
    def acoeff(self):
        """Computes the x values in Ax=b.
        """
        return np.dot(self.mmi,np.dot(self.mmatrix.T,self.bvec))
    #
    #
    #
    @property
    def chi2(self):
        """Computes the chi**2 value.
        """
        return np.sum((np.dot(self.mmatrix,self.acoeff) - self.bvec)**2)
    #
    #
    #
    @property
    def yfit(self):
        """Computes the best fit.
        """
        return np.dot(self.amatrix,self.acoeff)
    #
    #
    #
    @property
    def dof(self):
        """Computes the degrees of freedom.
        """
        return (self.sqivar > 0).sum() - self.nstar
    #
    #
    #
    @property
    def covar(self):
        """Computes the covariance matrix.
        """
        wwt = self.ww.copy()
        wwt[self.ww>0] = 1.0/self.ww[self.ww>0]
        covar = np.zeros((self.nstar,self.nstar),dtype=self.ww.dtype)
        for i in range(self.nstar):
            for j in range(i+1):
                covar[i,j] = np.sum(wwt * self.vv[:,i] * self.vv[:,j])
                covar[j,i] = covar[i,j]
        return covar
    #
    #
    #
    @property
    def var(self):
        """Computes the variances of the fit parameters.
        """
        return np.diag(self.covar)
