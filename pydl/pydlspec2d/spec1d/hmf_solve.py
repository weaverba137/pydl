# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from __future__ import print_function
#
import numpy as np
#
#
#
def model(a,g):
    """Compute the model."""
    return np.dot(a,g)
#
#
#
def resid(a,g,spectra):
    """Compute residuals."""
    return spectra - model(a,g)
#
#
#
def chi(a,g,spectra,invvar):
    """Compute chi, the scaled residual."""
    return resid(a,g,spectra)*np.sqrt(invvar)
#
#
#
def penalty(g,epsilon=None):
    """Compute penalty for non-smoothness."""
    if epsilon is None:
        return 0.0
    return epsilon*np.sum(np.diff(g)**2)
#
#
#
def badness(a,g,spectra,invvar,epsilon=None):
    """Compute chi**2."""
    return np.sum(chi(a,g,spectra,invvar)**2) + penalty(g,epsilon)
#
#
#
def normbase(g):
    """Apply standard component normalization."""
    return np.sqrt((g**2).mean(1))
#
#
#
def astep(spectra,invvar,g):
    """Update for coefficients at fixed component spectra."""
    from numpy.linalg import solve
    N,M = spectra.shape
    K,M = g.shape
    a = np.zeros((N,K),dtype=g.dtype)
    for i in range(N):
        Gi = np.zeros((K,K),dtype=g.dtype)
        for k in range(K):
            for kp in range(k,K):
                Gi[k,kp] = np.sum(g[k,:]*g[kp,:]*invvar[i,:])
                if kp > k:
                    Gi[kp,k] = Gi[k,kp]
        Fi = np.dot(g,spectra[i,:]*invvar[i,:])
        a[i,:] = solve(Gi,Fi)
    return a
#
#
#
def gstep(oldg,spectra,invvar,a,epsilon=None):
    """Update for component spectra at fixed coefficients."""
    from numpy.linalg import solve
    N,M = spectra.shape
    N,K = a.shape
    g = np.zeros((K,M),dtype=a.dtype)
    e = np.zeros(oldg.shape,dtype=oldg.dtype)
    d = np.zeros((K,K,M),dtype=a.dtype)
    if epsilon is not None and epsilon > 0:
        foo = epsilon*np.eye(K,dtype=a.dtype)
        for l in range(M):
            d[:,:,l] = foo
            if l > 0 and l < M-1:
                d[:,:,l] *= 2
        # d[:,:,0] = foo
        # d[:,:,1:M-1] = 2*foo
        # d[:,:,M-1] = foo
        e[:,0] = epsilon*oldg[:,1]
        e[:,1:M-1] = epsilon*(oldg[:,0:M-2] + oldg[:,2:M])
        e[:,M-1] = epsilon*oldg[:,M-2]
    for j in range(M):
        Aj = np.zeros((K,K),dtype=a.dtype)
        for k in range(K):
            for kp in range(k,K):
                Aj[k,kp] = np.sum(a[:,k]*a[:,kp]*invvar[:,j])
                if kp > k:
                    Aj[kp,k] = Aj[k,kp]
        Aj += d[:,:,j]
        Fj = np.dot(a.T,spectra[:,j]*invvar[:,j]) + e[:,j]
        g[:,j] = solve(Aj,Fj)
    return g
#
#
#
def astepnn(a,spectra,invvar,g):
    """Non-negative update for coefficients at fixed component spectra."""
    numerator = np.dot(spectra*invvar,g.T)
    denominator = np.dot(np.dot(a,g)*invvar,g.T)
    return a*(numerator/denominator)
#
#
#
def gstepnn(g,spectra,invvar,a,epsilon=None):
    """Non-negative update for component spectra at fixed coefficients."""
    K,M = g.shape
    numerator = np.dot(a.T,(spectra*invvar))
    if epsilon is not None and epsilon > 0:
        e = np.zeros(g.shape,dtype=g.dtype)
        e[:,0] = epsilon*g[:,1]
        e[:,1:M-1] = epsilon*(g[:,0:M-2] + g[:,2:M])
        e[:,M-1] = epsilon*g[:,M-2]
        numerator += e
    denominator = np.dot(a.T,np.dot(a,g)*invvar)
    if epsilon is not None and epsilon > 0:
        d = epsilon*g.copy()
        d[:,1:M-1] *= 2
        denominator += d
    return g*(numerator/denominator)
#
#
#
def reorder(a,g):
    """Reorder and rotate basis analogous to PCA."""
    from numpy.linalg import eigh
    l,U = eigh(np.dot(a.T,a))
    return (np.dot(a,U),np.dot(U.T,g))
#
#
#
def hmf_solve(spectra,invvar,K=4,nonnegative=False,epsilon=None):
    """Handle the HMF iteration, analogous to pca_solve.

    Parameters
    ----------
    spectra : numpy.array
        The training spectra.
    invvar : numpy.array
        The inverse variance for each pixel in the spectra.
    K : int, optional
        The number of dimensions of the factorization (default 4).
    nonnegative : bool, optional
        Set this to ``True`` to perform nonnegative HMF.
    epsilon : float, optional
        Regularization parameter.  Set to any non-zero float value to turn it on.

    Returns
    -------
    a,g : tuple of numpy.array
        The fitting coefficients and fitted functions, respectively.
    """
    from scipy.cluster.vq import kmeans, whiten
    N,M = spectra.shape
    #
    # Make spectra non-negative
    #
    if nonnegative:
        spectra[spectra<0] = 0
        invvar[spectra<0] = 0
    #
    # Detect very bad columns
    #
    si = spectra*invvar
    if (spectra.sum(0) == 0).any():
        raise ValueError("Columns of zeros detected in spectra!")
    if (invvar.sum(0) == 0).any():
        raise ValueError("Columns of zeros detected in invvar!")
    if (si.sum(0) == 0).any():
        raise ValueError("Columns of zeros detected in spectra*invvar!")
    #
    # Initialize g matrix with kmeans
    #
    whitespectra = whiten(spectra)
    g,foo = kmeans(whitespectra,K)
    g /= np.repeat(normbase(g),M).reshape(g.shape)
    #
    # Initialize a matrix
    #
    a = np.outer(np.sqrt((spectra**2).mean(1)),np.repeat(1.0/K,K))
    if nonnegative:
        for k in range(128):
            a = astepnn(a,spectra,invvar,g)
    #
    # Number of iterations.
    #
    if nonnegative:
        n_iter = 2048
    else:
        n_iter = 16
    #
    # Iterate!
    #
    for m in range(n_iter):
        print(m)
        if nonnegative:
            a = astepnn(a,spectra,invvar,g)
            g = gstepnn(g,spectra,invvar,a,epsilon)
        else:
            a = astep(spectra,invvar,g)
            g = gstep(g,spectra,invvar,a,epsilon)
            a,g = reorder(a,g)
        norm = normbase(g)
        g /= np.repeat(norm,M).reshape(g.shape)
        a = (a.T*np.repeat(norm,N).reshape(K,N)).T
        # print(badness(a,g,spectra,invvar,epsilon))
    return (a,g)
