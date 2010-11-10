#
#
#
def number_of_fibers(plate,**kwargs):
    """Returns the total number of fibers per plate.
    """
    import os
    import os.path
    import pyfits
    from pydlspec2d.spec1d import latest_mjd
    mjd = latest_mjd(plate,**kwargs)
    nfiber = np.zeros(mjd.size,dtype='i4')
    #
    # Not all BOSS plates have 1000 fibers
    #
    if 'path' in kwargs:
        platelistpath = os.path.join(kwargs['path'],'platelist.fits')
    else:
        platelistpath = os.path.join(os.getenv('BOSS_SPECTRO_REDUX'),'platelist.fits')
    platelist = pyfits.open(platelistpath)
    platentotal = platelist[1].data.field('N_TOTAL')
    plateplate = platelist[1].data.field('PLATE')
    platemjd = platelist[1].data.field('MJD')
    platerun2d = platelist[1].data.field('RUN2D')
    platerun1d = platelist[1].data.field('RUN1D')
    platelist.close()
    if 'run2d' in kwargs:
        run2d = kwargs['run2d']
    else:
        run2d = os.getenv('RUN2D')
    if 'run1d' in kwargs:
        run1d = kwargs['run1d']
    else:
        run1d = os.getenv('RUN1D')
    for k in range(mjd.size):
        nfiber[k] = platentotal[(plateplate==plate[k]) & (platemjd==mjd[k]) &
        (platerun2d==run2d) & (platerun1d==run1d)]
    #
    # SDSS-I,II plates
    #
    nfiber[mjd < 55025] = 640
    return nfiber
