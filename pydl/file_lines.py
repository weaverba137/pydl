# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-


def file_lines(path, compress=False):
    """Replicates the IDL ``FILE_LINES()`` function.

    Given a path to a file name or a list of such paths, returns the number of
    lines in the file(s).

    Parameters
    ----------
    path : :class:`str` or :class:`list` of :class:`str`
        Path to a file.  Can be a list of paths.
    compress : :class:`bool`, optional
        If set to ``True``, assumes that all files in `path` are GZIP
        compressed.

    Returns
    -------
    :class:`int` or :class:`list` of :class:`int`
        The number of lines in `path`.  Returns a list of lengths if a list of
        files is supplied.

    Notes
    -----
    The ``/NOEXPAND_PATH`` option in IDL's ``FILE_LINES()`` is not implemented.

    References
    ----------
    http://www.harrisgeospatial.com/docs/file_lines.html

    Examples
    --------
    >>> from pydl import file_lines
    >>> from os.path import dirname, join
    >>> file_lines(join(dirname(__file__),'tests','t','this-file-contains-42-lines.txt'))
    42
    """
    scalar = False
    if isinstance(path, (str,)):
        working_path = [path]
        scalar = True
    else:
        working_path = path
    lines = list()
    for filename in working_path:
        if compress:
            #
            # gzip in Python 2.6 can't use a context manager.
            #
            import gzip
            f = gzip.open(filename)
            lines.append(len(f.readlines()))
            f.close()
        else:
            with open(filename) as f:
                lines.append(len(f.readlines()))
    if scalar:
        return lines[0]
    else:
        return lines
