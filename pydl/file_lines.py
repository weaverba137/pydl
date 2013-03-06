# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def file_lines(path,compress=False):
    """Replicates the IDL FILE_LINES() function.

    Given a path to a file name or a list of such paths, returns the number of
    lines in the file(s).

    Parameters
    ----------
    path : str or list of str
        Path to a file.  Can be a list of paths.
    compress : bool, optional
        If set to ``True``, assumes that all files in `path` are GZIP compressed.

    Returns
    -------
    file_lines : int or list of int
        The number of lines in `path`.  Returns a list of lengths if a list of
        files is supplied.

    Notes
    -----
    The ``/NOEXPAND_PATH`` option in IDL's FILE_LINES() is not implemented.

    References
    ----------
    http://www.exelisvis.com/docs/FILE_LINES.html

    Examples
    --------
    >>> pydl.file_lines('my-one-line-file.txt')
    1
    """
    scalar = False
    if isinstance(path, str):
        working_path = [path]
        scalar = True
    else:
        working_path = path
    lines = list()
    for filename in working_path:
        if compress:
            import gzip
            with gzip.open(filename) as f:
                lines.append(len(f.readlines()))
        else:
            with open(filename) as f:
                lines.append(len(f.readlines()))
    if scalar:
        return lines[0]
    else:
        return lines
