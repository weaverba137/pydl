.. _pydl.goddard:

==================================
Goddard Utilities (`pydl.goddard`)
==================================

Introduction
++++++++++++

This package provides functionality similar to the `The IDL® Astronomy User's Libary`_,
sometimes called the "Goddard Utilities", maintained by Wayne Landsman and
distributed with idlutils_.

The Goddard package is itself divided into a number of subpackages.  Below we list
the subpackages and the usability of the PyDL equivalent.
The readiness levels are defined as:

Not Applicable (NA)
    Another built-in or numpy/scipy/astropy package completely replaces this.  No point in implementing.
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
astro         Rudimentary     Foo
astrom        None            Foo
coyote        None            Foo
database      None            Foo
disk_io       None            Foo
fits          NA              Use :mod:`astropy.io.fits`.
fits_bintable NA              Use :mod:`astropy.io.fits`.
fits_table    NA              Use :mod:`astropy.io.fits`.
idlphot       None            Foo
image         None            Foo
jhuapl        None            Foo
markwardt     None            Foo
math          Rudimentary     Foo
misc          Rudimentary     Foo
plot          None            Foo
robust        None            Foo
sdas          None            Foo
sockets       None            Foo
structure     None            Foo
tv            None            Foo
============= =============== ===================================================

.. _`The IDL® Astronomy User's Libary`: http://idlastro.gsfc.nasa.gov/
.. _idlutils: http://www.sdss.org/dr14/software/idlutils/

API
+++

.. automodapi:: pydl.goddard

.. automodapi:: pydl.goddard.astro

.. automodapi:: pydl.goddard.math

.. automodapi:: pydl.goddard.misc
