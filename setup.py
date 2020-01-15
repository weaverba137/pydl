#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Needed to import local helpers when pyproject.toml exists.
import sys
sys.path.append(".")

# Ensure that astropy-helpers is available
import ah_bootstrap  # noqa

from setuptools import setup

from astropy_helpers.setup_helpers import register_commands
# from distutils.command.sdist import sdist as DistutilsSdist

# Create a dictionary with setup command overrides. Note that this gets
# information about the package (name and version) from the setup.cfg file.
cmdclass = register_commands()
# cmdclass = {'sdist': DistutilsSdist}
cmdclass.pop('build_sphinx', None)
cmdclass.pop('build_docs', None)
# cmdclass.pop('sdist', None)
cmdclass.pop('build_ext', None)
# print(cmdclass)

setup(use_scm_version={"write_to": "pydl/version.py"}, cmdclass=cmdclass)
