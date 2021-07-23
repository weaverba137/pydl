====
PyDL
====

Introduction
++++++++++++


This package consists of Python_ replacements for functions that are part of
the `IDL®`_ built-in library or part of astronomical `IDL®`_ libraries.
The emphasis is on reproducing results of the astronomical library functions.
Only the bare minimum of `IDL®`_ built-in functions are implemented to support this.

There are four astronomical libraries targeted:

* idlutils_ : a general suite of tools heavily used by SDSS_.
* `Goddard utilities`_ : The `IDL®`_ Astronomy User's Libary, maintained by Wayne Landsman and distributed with idlutils_.
* idlspec2d_ : tools for working with SDSS_, BOSS_ and eBOSS_ spectroscopic data.
* photoop_ : tools for working with SDSS_ imaging data.

This package affiliated with the astropy_ project and is registered with PyPI_.

IDL® is a registered trademark of `Harris Geospatial Solutions`_.

Components
++++++++++

Most of the functionality of PyDL is in sub-packages.

.. toctree::
   :maxdepth: 1

   pydlutils.rst
   goddard.rst
   pydlspec2d.rst
   photoop.rst

Other Notes
+++++++++++

.. toctree::
   :maxdepth: 1

   templates.rst
   window.rst
   changes.rst
   todo.rst
   credits.rst
   licenses.rst

Base API
++++++++

.. automodapi:: pydl
   :no-inheritance-diagram:


.. _Python: https://www.python.org
.. _`IDL®`: https://www.l3harrisgeospatial.com/Software-Technology/IDL
.. _idlutils: https://www.sdss.org/dr16/software/idlutils/
.. _SDSS: https://www.sdss.org
.. _`Goddard utilities`: https://idlastro.gsfc.nasa.gov/
.. _idlspec2d: https://svn.sdss.org/public/repo/eboss/idlspec2d/trunk/
.. _BOSS: https://www.sdss.org/surveys/boss/
.. _eBOSS: https://www.sdss.org/surveys/eboss/
.. _photoop: https://svn.sdss.org/public/repo/sdss/photoop/trunk/
.. .. _astropy: http://www.astropy.org
.. _PyPI: https://pypi.org/project/pydl/
.. _`PyDL on Read the Docs`: https://pydl.readthedocs.io/en/latest/
.. _SDSS-III: http://www.sdss3.org
.. _`svn repository`: https://www.sdss.org/dr16/software/products/
.. _GitHub: https://github.com
.. _`Harris Geospatial Solutions`: https://www.l3harrisgeospatial.com/
