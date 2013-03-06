"""
====
pydl
====

This docstring is supplied to the long description in setup.py.
"""
from file_lines import file_lines
from pcomp import pcomp
from smooth import smooth
from uniq import uniq

class PydlException(Exception):
    pass

__all__ = ['PydlException']
