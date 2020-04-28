.. _pydl.goddard:

==================================
Goddard Utilities (`pydl.goddard`)
==================================

Introduction
++++++++++++

This package provides functionality similar to the `The IDL® Astronomy User's Libary`_,
sometimes called the "Goddard Utilities", maintained by Wayne Landsman and
distributed with idlutils_.

In general, functions that are needed by :mod:`pydl.pydlutils` or :mod:`pydl.pydlspec2d`
are implemented, while functions that are *not* needed have much lower priority.

The Goddard package is itself divided into a number of subpackages.  Below we list
the subpackages and the usability of the PyDL equivalent.
The readiness levels are defined as:

Obsolete
    No point in implementing because the purpose of the code lapsed many years ago.
Not Applicable (NA)
    No point in implementing because another built-in or numpy/scipy/astropy package completely replaces this.
None
    Not (yet) implemented at all.
Rudimentary
    Only a few functions are implemented.
Fair
    Enough functions are implemented to be useful, but some are missing.
Good
    Pretty much anything you could do with the Goddard code you can do with the equivalent here.

============= =============== ===================================================
Subpackage    Readiness Level Comments
============= =============== ===================================================
astro         Rudimentary     General astronomical utility functions.
astrom        NA              Tools for manipulating WCS data in FITS headers.  Use :mod:`astropy.io.fits` and :mod:`astropy.wcs`.
coyote        NA              The `Coyote library`_ for plotting and graphics developed by David Fanning.
database      None            Allows access to IDL-specific databases.
disk_io       None            Provides access to IRAF image (``.imh``) files and AJ/ApJ-style tables.
fits          NA              Use :mod:`astropy.io.fits`.
fits_bintable NA              Use :mod:`astropy.io.fits`.
fits_table    NA              Use :mod:`astropy.io.fits`.
idlphot       None            Adapted from an early version of DAOPHOT.
image         None            Generic image processing functions, including convolution/deconvolution.
jhuapl        None            `Functions from the JHU Applied Physics Lab`_.
markwardt     None            Levenberg-Marquardt least-squares minimization.
math          Rudimentary     Generic mathematical functions.  Many are implemented in numpy or scipy.
misc          Rudimentary     General utility functions that do not involve astronomy specifically.
plot          NA              Functions that supplement the built-in IDL plotting capabilities.
robust        None            Robust statistical fitting procedures.
sdas          None            Provides access to `STDAS/GEIS`_ image files.
sockets       NA              Functions for performing web queries in IDL.  Use astroquery_.
structure     NA              Tools for manipulating IDL data structures. Use :class:`numpy.recarray`.
tv            NA              Functions for manipulating IDL image displays.
============= =============== ===================================================

.. _`The IDL® Astronomy User's Libary`: https://idlastro.gsfc.nasa.gov/
.. _idlutils: https://www.sdss.org/dr16/software/idlutils/
.. _`Coyote library`: http://www.idlcoyote.com/
.. _`Functions from the JHU Applied Physics Lab`: http://fermi.jhuapl.edu/s1r/idl/idl.html
.. _`STDAS/GEIS`: http://www.stsci.edu/instruments/wfpc2/Wfpc2_dhb/intro_ch24.html#1905747
.. _astroquery: https://astroquery.readthedocs.io/en/latest/

API
+++

.. automodapi:: pydl.goddard

.. automodapi:: pydl.goddard.astro
    :skip: time

.. automodapi:: pydl.goddard.math
    :skip: legendre

.. automodapi:: pydl.goddard.misc
