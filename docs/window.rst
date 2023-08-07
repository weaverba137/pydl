====================================
The SDSS Photometric Window Function
====================================

The "window" function is the function that answers the question,
"Is this point (RA, Dec) in an SDSS imaging survey image?"  In practice,
it's easiest to think of the window function as a big look-up table.
The determination of the window function is intrinsically intwined with
the "resolve" algorithm, which determines the highest quality ("primary")
SDSS image that covers a region of the sky.  For further low-level details,
see the `SDSS documentation`_ of the resolve algorithm.

.. _`SDSS documentation`: https://www.sdss4.org/dr16/algorithms/resolve/

In the PyDL package, the :mod:`~pydl.photoop.window` and :mod:`~pydl.pydlutils.mangle`
modules are used to access the window function.  In this document we will
concentrate on these questions:

* Is a given point or set of points in the SDSS imaging region?
* If a point *is* in the SDSS imaging, which specific image contains it?

.. All you really need to do is read window_unified.fits with read_fits_polygons.
