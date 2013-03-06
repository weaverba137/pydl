#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import sys
import imp
try:
    # This incantation forces distribute to be used (over setuptools) if it is
    # available on the path; otherwise distribute will be downloaded.
    import pkg_resources
    distribute = pkg_resources.get_distribution('distribute')
    if pkg_resources.get_distribution('setuptools') != distribute:
        sys.path.insert(1, distribute.location)
        distribute.activate()
        imp.reload(pkg_resources)
except:  # There are several types of exceptions that can occur here
    from distribute_setup import use_setuptools
    use_setuptools()

import glob
import os
from setuptools import setup, find_packages

#A dirty hack to get around some early import/configurations ambiguities
if sys.version_info[0] >= 3:
    import builtins
else:
    import __builtin__ as builtins
builtins._PACKAGE_SETUP_ = True

import astropy
from astropy.setup_helpers import (register_commands, adjust_compiler,
                                   filter_packages, update_package_files,
                                   get_debug_option)
from astropy.version_helpers import get_git_devstr, generate_version_py

# Use this import to get the docstring
import pydl

# Set affiliated package-specific settings
PACKAGENAME = 'pydl'
DESCRIPTION = 'Astropy affiliated package'
LONG_DESCRIPTION = pydl.__doc__
AUTHOR = 'Benjamin Alan Weaver'
AUTHOR_EMAIL = 'benjamin.weaver@nyu.edu'
LICENSE = 'BSD'
URL = 'http://github.com/weaverba137/pydl'

# VERSION should be PEP386 compatible (http://www.python.org/dev/peps/pep-0386)
VERSION = '0.2.dev'

# Indicates if this version is a release version
RELEASE = 'dev' not in VERSION

if not RELEASE:
    VERSION += get_git_devstr(False)

# Populate the dict of setup command overrides; this should be done before
# invoking any other functionality from distutils since it can potentially
# modify distutils' behavior.
cmdclassd = register_commands(PACKAGENAME, VERSION, RELEASE)

# Adjust the compiler in case the default on this platform is to use a
# broken one.
adjust_compiler(PACKAGENAME)

# Freeze build information in version.py
generate_version_py(PACKAGENAME, VERSION, RELEASE, get_debug_option())

# Use the find_packages tool to locate all packages and modules
packagenames = filter_packages(find_packages())

# Treat everything in scripts except README.rst as a script to be installed
scripts = [fname for fname in glob.glob(os.path.join('scripts', '*'))
           if fname != 'README.rst']

# Additional C extensions that are not Cython-based should be added here.
extensions = []

# A dictionary to keep track of all package data to install
package_data = {PACKAGENAME: ['data/*']}

# A dictionary to keep track of extra packagedir mappings
package_dirs = {}

# Update extensions, package_data, packagenames and package_dirs from
# any sub-packages that define their own extension modules and package
# data.  See the docstring for setup_helpers.update_package_files for
# more details.
update_package_files(PACKAGENAME, extensions, package_data, packagenames,
                     package_dirs)


setup(name=PACKAGENAME,
      version=VERSION,
      description=DESCRIPTION,
      packages=packagenames,
      package_data=package_data,
      package_dir=package_dirs,
      ext_modules=extensions,
      scripts=scripts,
      requires=['astropy'],
      install_requires=['astropy'],
      provides=[PACKAGENAME],
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      license=LICENSE,
      url=URL,
      long_description=LONG_DESCRIPTION,
      cmdclass=cmdclassd,
      zip_safe=False,
      use_2to3=True
      )

#~ setup(
    #~ name='pydl',
    #~ version='0.1',
    #~ author='Benjamin Alan Weaver',
    #~ author_email='benjamin.weaver@nyu.edu',
    #~ packages=[
        #~ 'pydl',
        #~ 'pydl.test',
        #~ 'pydlutils',
        #~ 'pydlutils.bspline',
        #~ 'pydlutils.goddard',
        #~ 'pydlutils.goddard.astro',
        #~ 'pydlutils.goddard.fits',
        #~ 'pydlutils.goddard.math',
        #~ 'pydlutils.image',
        #~ 'pydlutils.mangle',
        #~ 'pydlutils.math',
        #~ 'pydlutils.misc',
        #~ 'pydlutils.sdss',
        #~ 'pydlutils.test',
        #~ 'pydlspec2d',
        #~ 'pydlspec2d.spec1d',
        #~ 'pydlspec2d.spec2d',
        #~ ],
    #~ url='http://pypi.python.org/pypi/pydl/',
    #~ license='LICENSE.txt',
    #~ description='Replacement for IDL',
    #~ long_description=open('README.txt').read(),
    #~ classifiers = [ 'Development Status :: 2 - Pre-Alpha',
        #~ 'Environment :: Console',
        #~ 'Intended Audience :: Science/Research',
        #~ 'License :: OSI Approved :: GNU General Public License (GPL)',
        #~ 'Operating System :: POSIX :: Linux',
        #~ 'Programming Language :: Python',
        #~ 'Topic :: Scientific/Engineering',
        #~ 'Topic :: Scientific/Engineering :: Astronomy',
        #~ ],
#~ )
#~
