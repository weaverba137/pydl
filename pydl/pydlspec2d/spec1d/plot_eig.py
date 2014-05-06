# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from __future__ import print_function
#
def plot_eig(filename,title='Unknown'):
    """Plot spectra from an eigenspectra/template file.

    Parameters
    ----------
    filename : str
        Name of a FITS file containing eigenspectra/templates.
    title : str, optional
        Title to put on the plot.

    Returns
    -------
    None
    """
    from astropy.io import fits as pyfits
    import pylab
    #
    # Set title based on filename
    #
    if title == 'Unknown':
        if filename.find('Gal') > 0:
            title = 'Galaxies: Eigenspectra'
        elif filename.find('QSO') > 0:
            title = 'QSOs: Eigenspectra'
        elif filename.find('Star') > 0:
            title = 'Stars: Eigenspectra'
        elif filename.find('CVstar') > 0:
            title = 'CV Stars: Eigenspectra'
        else:
            print('Unknown template type!')
    base,ext = filename.split('.')
    spectrum = pyfits.open(filename)
    newloglam0 = spectrum[0].header['COEFF0']
    objdloglam = spectrum[0].header['COEFF1']
    spectro_data = spectrum[0].data
    spectrum.close()
    (neig, ndata) = spectro_data.shape
    newloglam = pylab.arange(ndata) * objdloglam + newloglam0
    lam = 10.0**newloglam
    fig = pylab.figure(dpi=100)
    ax = fig.add_subplot(111)
    colorvec = ['k','r','g','b','m','c']
    for l in range(neig):
        p = ax.plot(lam,spectro_data[l,:],'%s-'%colorvec[l%len(colorvec)],linewidth=1)
    ax.set_xlabel(r'Wavelength [$\AA$]')
    ax.set_ylabel('Flux [Arbitrary Units]')
    ax.set_title(title)
    # ax.set_xlim([3500.0,10000.0])
    # ax.set_ylim([-400.0,500.0])
    # fig.savefig(base+'.zoom.png')
    fig.savefig(base+'.png')
    pylab.close(fig)
    return
