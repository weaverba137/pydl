# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""This module corresponds to the math directory in idlutils.
"""
import numpy as np
import astropy.utils as au


class computechi2(object):
    """Solve the linear set of equations :math:`A x = b` using SVD.

    The attributes of this class are all read-only properties, implemented
    with :class:`~astropy.utils.decorators.lazyproperty`.

    Parameters
    ----------
    bvec : :class:`numpy.ndarray`
        The :math:`b` vector in :math:`A x = b`. This vector has length
        :math:`N`.
    sqivar : :class:`numpy.ndarray`
        The reciprocal of the errors in `bvec`.  The name comes from the square
        root of the inverse variance, which is what this is.
    amatrix : :class:`numpy.ndarray`
        The matrix :math:`A` in :math:`A x = b`.
        The shape of this matrix is (:math:`N`, :math:`M`).
    """

    def __init__(self, bvec, sqivar, amatrix):
        """Initialize the object and perform initial computations.
        """
        from numpy.linalg import svd
        #
        # Save the inputs
        #
        # self.bvec = bvec
        self.sqivar = sqivar
        self.amatrix = amatrix
        if len(amatrix.shape) > 1:
            self.nstar = amatrix.shape[1]
        else:
            self.nstar = 1
        self.bvec = bvec * sqivar
        self.mmatrix = self.amatrix * np.tile(sqivar, self.nstar).reshape(
                       self.nstar, bvec.size).transpose()
        mm = np.dot(self.mmatrix.T, self.mmatrix)
        self.uu, self.ww, self.vv = svd(mm, full_matrices=False)
        self.mmi = np.dot((self.vv.T / np.tile(self.ww, self.nstar).reshape(
                   self.nstar, self.nstar)), self.uu.T)
        return

    @au.lazyproperty
    def acoeff(self):
        """(:class:`~numpy.ndarray`) The fit parameters, :math:`x`,
        in :math:`A x = b`. This vector has length :math:`M`.
        """
        return np.dot(self.mmi, np.dot(self.mmatrix.T, self.bvec))

    @au.lazyproperty
    def chi2(self):
        """(:class:`float <numpy.generic>`) The :math:`\chi^2` value of the fit.
        """
        return np.sum((np.dot(self.mmatrix, self.acoeff) - self.bvec)**2)

    @au.lazyproperty
    def yfit(self):
        """(:class:`~numpy.ndarray`) The evaluated best-fit at each point.
        This vector has length :math:`N`.
        """
        return np.dot(self.amatrix, self.acoeff)

    @au.lazyproperty
    def dof(self):
        """(:class:`int <numpy.generic>`) The degrees of freedom of the fit.
        This is the number of values of `bvec` that have `sqivar` > 0 minus
        the number of fit paramaters, which is equal to :math:`M`.
        """
        return (self.sqivar > 0).sum() - self.nstar

    @au.lazyproperty
    def covar(self):
        """(:class:`~numpy.ndarray`) The covariance matrix.
        The shape of this matrix is (:math:`M`, :math:`M`).
        """
        wwt = self.ww.copy()
        wwt[self.ww > 0] = 1.0/self.ww[self.ww > 0]
        covar = np.zeros((self.nstar, self.nstar), dtype=self.ww.dtype)
        for i in range(self.nstar):
            for j in range(i + 1):
                covar[i, j] = np.sum(wwt * self.vv[:, i] * self.vv[:, j])
                covar[j, i] = covar[i, j]
        return covar

    @au.lazyproperty
    def var(self):
        """(:class:`~numpy.ndarray`) The variances of the fit.
        This is identical to the diagonal of the covariance matrix.
        This vector has length :math:`M`.
        """
        return np.diag(self.covar)


def djs_median(array, dimension=None, width=None, boundary='none'):
    """Compute the median of an array.

    Use a filtering box or collapse the image along one dimension.

    Parameters
    ----------
    array : :class:`numpy.ndarray`
        input array
    dimension : :class:`int`, optional
        Compute the median over this dimension. It is an error to specify both
        `dimension` and `width`.
    width : :class:`int`, optional
        Width of the median window. In general, this should be an odd
        integer.  It is an error to specify both `dimension` and `width`.
    boundary : { 'none', 'reflect', 'nearest', 'wrap' }, optional
        Boundary condition to impose.  'none' means no filtering is done within
        `width`/2 of the boundary.  'reflect' means reflect pixel values around the
        boundary. 'nearest' means use the values of the nearest boundary pixel.
        'wrap' means wrap pixel values around the boundary. 'nearest' and 'wrap'
        are not implemented.

    Returns
    -------
    :class:`numpy.ndarray`
        The output.  If neither `dimension` nor `width` are set, this is a scalar
        value, just the output of ``numpy.median()``.  If `dimension` is set,
        then the result simply ``numpy.median(array,dimension)``.
        If `width` is set, the result has the same shape as the input array.
    """
    from ..median import median
    if dimension is None and width is None:
        return np.median(array)
    elif width is None:
        return np.median(array, axis=dimension)
    elif dimension is None:
        if width == 1:
            return array
        if boundary == 'none':
            if array.ndim == 1:
                return median(array, width)
            elif array.ndim == 2:
                return median(array, width)
            else:
                raise ValueError('Unsupported number of dimensions with ' +
                                'this boundary condition.')
        elif boundary == 'reflect':
            padsize = int(np.ceil(width/2.0))
            if array.ndim == 1:
                bigarr = np.zeros(array.shape[0]+2*padsize, dtype=array.dtype)
                bigarr[padsize:padsize+array.shape[0]] = array
                bigarr[0:padsize] = array[0:padsize][::-1]
                bigarr[padsize+array.shape[0]:padsize*2+array.shape[0]] = (
                        array[array.shape[0]-padsize:array.shape[0]][::-1])
                f = median(bigarr, width)
                medarray = f[padsize:padsize+array.shape[0]]
                return medarray
            elif array.ndim == 2:
                bigarr = np.zeros((array.shape[0]+2*padsize,
                                    array.shape[1]+2*padsize),
                                    dtype=array.dtype)
                bigarr[padsize:padsize+array.shape[0], padsize:padsize+array.shape[1]] = array
                # Copy into top + bottom
                bigarr[0:padsize, padsize:array.shape[1]+padsize] = array[0:padsize, :][::-1, :]
                bigarr[array.shape[0]+padsize:bigarr.shape[0], padsize:array.shape[1]+padsize] = array[array.shape[0]-padsize:array.shape[0], :][::-1, :]
                # Copy into left + right
                bigarr[padsize:array.shape[0]+padsize, 0:padsize] = array[:, 0:padsize][:, ::-1]
                bigarr[padsize:array.shape[0]+padsize, array.shape[1]+padsize:bigarr.shape[1]] = array[:, array.shape[1]-padsize:array.shape[1]][:, ::-1]
                # Copy into top left
                bigarr[0:padsize, 0:padsize] = array[0:padsize, 0:padsize][::-1, ::-1]
                # Copy into top right
                bigarr[0:padsize, bigarr.shape[1]-padsize:bigarr.shape[1]] = array[0:padsize, array.shape[1]-padsize:array.shape[1]][::-1, ::-1]
                # Copy into bottom left
                bigarr[bigarr.shape[0]-padsize:bigarr.shape[0], 0:padsize] = array[array.shape[0]-padsize:array.shape[0], 0:padsize][::-1, ::-1]
                # Copy into bottom right
                bigarr[bigarr.shape[0]-padsize:bigarr.shape[0], bigarr.shape[1]-padsize:bigarr.shape[1]] = array[array.shape[0]-padsize:array.shape[0], array.shape[1]-padsize:array.shape[1]][::-1, ::-1]
                f = median(bigarr, min(width, array.size))
                medarray = f[padsize:array.shape[0]+padsize, padsize:array.shape[1]+padsize]
                return medarray
            else:
                raise ValueError('Unsupported number of dimensions with ' +
                                'this boundary condition.')
        elif boundary == 'nearest':
            raise ValueError('This boundary condition not implemented')
        elif boundary == 'wrap':
            raise ValueError('This boundary condition not implemented')
        else:
            raise ValueError('Unknown boundary condition.')
    else:
        raise ValueError('Invalid to specify both dimension & width.')


def djs_reject(data, model, outmask=None, inmask=None, sigma=None,
               invvar=None, lower=None, upper=None, maxdev=None,
               maxrej=None, groupdim=None, groupsize=None, groupbadpix=False,
               grow=0, sticky=False):
    """Routine to reject points when doing an iterative fit to data.

    Parameters
    ----------
    data : :class:`numpy.ndarray`
        The data
    model : :class:`numpy.ndarray`
        The model, must have the same number of dimensions as `data`.
    outmask : :class:`numpy.ndarray`, optional
        Output mask, generated by a previous call to `djs_reject`.  If not supplied,
        this mask will be initialized to a mask that masks nothing.  Although
        this parameter is technically optional, it will almost always be set.
    inmask : :class:`numpy.ndarray`, optional
        Input mask.  Bad points are marked with a value that evaluates to ``False``.
        Must have the same number of dimensions as `data`.
    sigma : :class:`numpy.ndarray`, optional
        Standard deviation of the data, used to reject points based on the values
        of `upper` and `lower`.
    invvar : :class:`numpy.ndarray`, optional
        Inverse variance of the data, used to reject points based on the values
        of `upper` and `lower`.  If both `sigma` and `invvar` are set, `invvar`
        will be ignored.
    lower : :class:`int` or :class:`float`, optional
        If set, reject points with data < model - lower * sigma.
    upper : :class:`int` or :class:`float`, optional
        If set, reject points with data > model + upper * sigma.
    maxdev : :class:`int` or :class:`float`, optional
        If set, reject points with abs(data-model) > maxdev.  It is permitted to
        set all three of `lower`, `upper` and `maxdev`.
    maxrej : :class:`int` or :class:`numpy.ndarray`, optional
        Maximum number of points to reject in this iteration.  If `groupsize` or
        `groupdim` are set to arrays, this should be an array as well.
    groupdim
        To be documented.
    groupsize
        To be documented.
    groupbadpix : :class:`bool`, optional
        If set to ``True``, consecutive sets of bad pixels are considered groups,
        overriding the values of `groupsize`.
    grow : :class:`int`, optional
        If set to a non-zero integer, N, the N nearest neighbors of rejected
        pixels will also be rejected.
    sticky : :class:`bool`, optional
        If set to ``True``, pixels rejected in one iteration remain rejected in
        subsequent iterations, even if the model changes.

    Returns
    -------
    :func:`tuple`
        A tuple containing a mask where rejected data values are ``False`` and
        a boolean value set to ``True`` if `djs_reject` believes there is no
        further rejection to be done.

    Raises
    ------
    ValueError
        If dimensions of various inputs do not match.
    """
    from .misc import djs_laxisnum
    #
    # Create outmask setting = 1 for good data.
    #
    if outmask is None:
        outmask = np.ones(data.shape, dtype='bool')
    else:
        if data.shape != outmask.shape:
            raise ValueError('Dimensions of data and outmask do not agree.')
    #
    # Check other inputs.
    #
    if model is None:
        if inmask is not None:
            outmask = inmask
        return (outmask, False)
    else:
        if data.shape != model.shape:
            raise ValueError('Dimensions of data and model do not agree.')
    if inmask is not None:
        if data.shape != inmask.shape:
            raise ValueError('Dimensions of data and inmask do not agree.')
    if maxrej is not None:
        if groupdim is not None:
            if len(maxrej) != len(groupdim):
                raise ValueError('maxrej and groupdim must have the same number of elements.')
        else:
            groupdim = []
        if groupsize is not None:
            if len(maxrej) != len(groupsize):
                raise ValueError('maxrej and groupsize must have the same number of elements.')
        else:
            groupsize = len(data)
    if sigma is None and invvar is None:
        if inmask is not None:
            igood = (inmask & outmask).nonzero()[0]
        else:
            igood = outmask.nonzero()[0]
        if len(igood > 1):
            sigma = np.std(data[igood] - model[igood])
        else:
            sigma = 0
    diff = data - model
    #
    # The working array is badness, which is set to zero for good points
    # (or points already rejected), and positive values for bad points.
    # The values determine just how bad a point is, either corresponding
    # to the number of sigma above or below the fit, or to the number
    # of multiples of maxdev away from the fit.
    #
    badness = np.zeros(outmask.shape, dtype=data.dtype)
    #
    # Decide how bad a point is according to lower.
    #
    if lower is not None:
        if sigma is not None:
            qbad = diff < (-lower * sigma)
            badness += ((-diff/(sigma + (sigma == 0))) > 0) * qbad
        else:
            qbad = (diff * np.sqrt(invvar)) < -lower
            badness += ((-diff * np.sqrt(invvar)) > 0) * qbad
    #
    # Decide how bad a point is according to upper.
    #
    if upper is not None:
        if sigma is not None:
            qbad = diff > (upper * sigma)
            badness += ((diff/(sigma + (sigma == 0))) > 0) * qbad
        else:
            qbad = (diff * np.sqrt(invvar)) > upper
            badness += ((diff * np.sqrt(invvar)) > 0) * qbad
    #
    # Decide how bad a point is according to maxdev.
    #
    if maxdev is not None:
        qbad = np.absolute(diff) > maxdev
        badness += np.absolute(diff) / maxdev * qbad
    #
    # Do not consider rejecting points that are already rejected by inmask.
    # Do not consider rejecting points that are already rejected by outmask,
    # if sticky is set.
    #
    if inmask is not None:
        badness *= inmask
    if sticky:
        badness *= outmask
    #
    # Reject a maximum of maxrej (additional) points in all the data, or
    # in each group as specified by groupsize, and optionally along each
    # dimension specified by groupdim.
    #
    if maxrej is not None:
        #
        # Loop over each dimension of groupdim or loop once if not set.
        #
        for iloop in range(max(len(groupdim), 1)):
            #
            # Assign an index number in this dimension to each data point.
            #
            if len(groupdim) > 0:
                yndim = len(ydata.shape)
                if groupdim[iloop] > yndim:
                    raise ValueError('groupdim is larger than the number of dimensions for ydata.')
                dimnum = djs_laxisnum(ydata.shape, iaxis=groupdim[iloop]-1)
            else:
                dimnum = 0
            #
            # Loop over each vector specified by groupdim. For example, if
            # this is a 2-D array with groupdim=1, then loop over each
            # column of the data.  If groupdim=2, then loop over each row.
            # If groupdim is not set, then use the whole image.
            #
            for ivec in range(max(dimnum)):
                #
                # At this point it is not possible that dimnum is not set.
                #
                indx = (dimnum == ivec).nonzero()[0]
                #
                # Within this group of points, break it down into groups
                # of points specified by groupsize, if set.
                #
                nin = len(indx)
                if groupbadpix:
                    goodtemp = badness == 0
                    groups_lower = (-1*np.diff(np.insert(goodtemp, 0, 1)) == 1).nonzero()[0]
                    groups_upper = (np.diff(np.append(goodtemp, 1)) == 1).nonzero()[0]
                    ngroups = len(groups_lower)
                else:
                    #
                    # The IDL version of this test makes no sense because
                    # groupsize will always be set.
                    #
                    if 'groupsize' in kwargs:
                        ngroups = nin/groupsize + 1
                        groups_lower = np.arange(ngroups, dtype='i4')*groupsize
                        foo = (np.arange(ngroups, dtype='i4')+1)*groupsize
                        groups_upper = np.where(foo < nin, foo, nin) - 1
                    else:
                        ngroups = 1
                        groups_lower = [0, ]
                        groups_upper = [nin - 1, ]
                for igroup in range(ngroups):
                    i1 = groups_lower[igroup]
                    i2 = groups_upper[igroup]
                    nii = i2 - i1 + 1
                    #
                    # Need the test that i1 != -1 below to prevent a crash
                    # condition, but why is it that we ever get groups
                    # without any points?  Because this is badly-written,
                    # that's why.
                    #
                    if nii > 0 and i1 != -1:
                        jj = indx[i1:i2+1]
                        #
                        # Test if too many points rejected in this group.
                        #
                        if np.sum(badness[jj] != 0) > maxrej[iloop]:
                            isort = badness[jj].argsort()
                            #
                            # Make the following points good again.
                            #
                            badness[jj[isort[0:nii-maxrej[iloop]]]] = 0
                        i1 += groupsize[iloop]
    #
    # Now modify outmask, rejecting points specified by inmask=0, outmask=0
    # if sticky is set, or badness > 0.
    #
    # print(badness)
    newmask = badness == 0
    # print(newmask)
    if grow > 0:
        rejects = newmask == 0
        if rejects.any():
            irejects = rejects.nonzero()[0]
            for k in range(1, grow):
                newmask[(irejects - k) > 0] = 0
                newmask[(irejects + k) < (data.shape[0]-1)] = 0
    if inmask is not None:
        newmask = newmask & inmask
    if sticky:
        newmask = newmask & outmask
    #
    # Set qdone if the input outmask is identical to the output outmask.
    #
    qdone = np.all(newmask == outmask)
    outmask = newmask
    return (outmask, qdone)


def find_contiguous(x):
    """Find the longest sequence of contiguous non-zero array elements.

    Parameters
    ----------
    x : :class:`numpy.ndarray`
        A 1d array. A dtype of bool is preferred although any dtype where the
        operation ``if x[k]:`` is well-defined should work.

    Returns
    -------
    :class:`list`
        A list of indices of the longest contiguous non-zero sequence.

    Examples
    --------
    >>> import numpy as np
    >>> from pydl.pydlutils.math import find_contiguous
    >>> find_contiguous(np.array([0,1,1,1,0,1,1,0,1]))
    [1, 2, 3]
    """
    contig = list()
    for k in range(x.size):
        if x[k]:
            if len(contig) == 0:
                contig.append([k])
            else:
                if k == contig[-1][-1]+1:
                    contig[-1].append(k)
                else:
                    contig.append([k])
    lengths = [len(c) for c in contig]
    longest = contig[lengths.index(max(lengths))]
    return longest
