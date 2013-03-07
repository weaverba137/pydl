#
#
#
def set_use_caps(x,cm,polygon_use_caps,**kwargs):
    """Set the bits in use_caps for a polygon.
    """
    import numpy as np
    from pydlutils.mangle import is_cap_used
    if 'add' in kwargs:
        use_caps = long(polygon_use_caps)
    else:
        use_caps = 0L
    if 'tol' not in kwargs:
        kwargs['tol'] = 1.0e-10
    if 'allow_neg_doubles' in kwargs:
        nd = kwargs['allow_neg_doubles']
    else:
        nd = False
    t2 = kwargs['tol']**2
    use_caps |= 2L**len(cm) - 1L
    if 'allow_doubles' not in kwargs:
        #
        # Check for doubles
        #
        for i in range(len(cm)):
            if is_cap_used(use_caps,i):
                for j in range(i+1,len(cm)):
                    if is_cap_used(use_caps,j):
                        if np.sum(x[i]-x[j])**2 < t2:
                            if ((np.absolute(cm[i]-cm[j]) < kwargs['tol']) or
                                ((cm[i] + cm[j]) < kwargs['tol'] and not nd)):
                                #
                                # Don't use
                                #
                                use_caps -= 2L**j
    return use_caps

