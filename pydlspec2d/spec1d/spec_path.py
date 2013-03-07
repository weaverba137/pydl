#
# $Id$
#
def spec_path(plate,**kwargs):
    """Return the path to spPlate files
    """
    import os
    import numpy as np
    if isinstance(plate,int) or isinstance(plate,long):
        platevec = np.zeros(1,dtype='i4') + plate
    else:
        platevec = plate
    if 'path' in kwargs:
        platepath = True
    else:
        platepath = False
        if 'topdir' in kwargs:
            topdir = kwargs['topdir']
        else:
            topdir = os.getenv('BOSS_SPECTRO_REDUX')
        if 'run2d' in kwargs:
            run2d = kwargs['run2d']
        else:
            run2d = os.getenv('RUN2D')
    paths = list()
    for p in platevec:
        if platepath:
            paths.append(kwargs['path'])
        else:
            paths.append("%s/%s/%04d" % (topdir,run2d,p))
    return paths

