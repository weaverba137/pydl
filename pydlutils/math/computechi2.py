def computechi2(y,sqivar,templates):
    """Solve the linear set of equations Ax=b using SVD.
    """
    import numpy as np
    if len(templates.shape) > 1:
        nstar = templates.shape[1]
    else:
        nstar = 1
    bvec = y*sqivar
    mmatrix = templates * np.tile(sqivar,nstar).reshape(nstar,y.size).transpose()
    mm = np.dot(mmatrix.T,mmatrix)
    uu,ww,vv = svd(mm,full_matrices=False)
    mmi = vv.T / np.tile(ww,nstar).reshape(nstar,nstar)
    mmi = np.dot(mmi,uu.T)
    r = dict()
    r['acoeff'] = np.dot(mmi,np.dot(mmatrix.T,bvec))
    r['chi2'] = np.sum((np.dot(mmatrix,r['acoeff']) - bvec)**2)
    r['yfit'] = np.dot(templates,r['acoeff'])
    r['dof'] =  (sqivar > 0).sum() - nstar
    wwt = ww.copy()
    wwt[ww>0] = 1.0/ww[ww>0]
    covar = np.zeros((nstar,nstar),dtype=ww.dtype)
    for i in range(nstar):
        for j in range(i+1):
            covar[i,j] = np.sum(wwt * vv[:,i] * vv[:,j])
            covar[j,i] = covar[i,j]
    r['covar'] = covar
    return r

