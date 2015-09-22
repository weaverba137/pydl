# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from __future__ import print_function
#
def struct_print(array,filename=None,formatcodes=None,alias=None,fdigit=5,ddigit=7,
    html=False,no_head=False,debug=False,silent=False):
    """Print a NumPy record array (analogous to an IDL structure) in a nice way.

    Parameters
    ----------
    array : :class:`numpy.ndarray`
        A record array to print.
    filename : :class:`str` or file-like, optional
        If supplied, write to this file.
    formatcodes : :class:`dict`, optional
        If supplied, use explicit format for certain columns.
    alias : :class:`dict`, optional
        If supplied, use this mapping of record array column names to printed column names.
    fdigit : :class:`int`, optional
        Width of 32-bit floating point columns, default 5.
    ddigit : :class:`int`, optional
        Width of 64-bit floating point columns, default 7.
    html : :class:`bool`, optional
        If ``True``, print an html table.
    no_head : :class:`bool`, optional
        If ``True``, *don't* print a header line.
    debug : :class:`bool`, optional
        If ``True``, print some extra debugging information.
    silent : :class:`bool`, optional
        If ``True``, do not print the table, just return it.
    Returns
    -------
    struct_print : :func:`tuple`
        A tuple containing a list of the lines in the table.  If `html` is ``True``,
        also returns a list of lines of CSS for formatting the table.

    Examples
    --------
    >>> import numpy as np
    >>> from pydl.pydlutils.misc import struct_print
    >>> struct_print(np.array([(1,2.34,'five'),(2,3.456,'seven'),(3,4.5678,'nine')],dtype=[('a','i4'),('bb','f4'),('ccc','S5')]),silent=True)
    (['a bb          ccc  ', '- ----------- -----', '1        2.34 five ', '2       3.456 seven', '3      4.5678 nine '], [])
    """
    import numpy as np
    from . import decode_mixed
    from .. import PydlutilsException
    f = None # This variable will store a file handle
    if filename is not None:
        if isinstance(filename,file):
            f = filename
        else:
            f = open(filename,'w')
    if html:
        headstart = '<tr><th>'
        headsep = '</th><th>'
        headend = '</th></tr>'
        colstart = '<tr><td>'
        colsep = '</td><td>'
        colend = '</td></tr>'
        css = [ '<style type="text/css">',
            'table {',
            '    border-collapse: collapse;',
            '}',
            'th {',
            '    padding: 2px;',
            '    text-align: right;',
            '    border: 1px solid black;',
            '    font-weight: bold;',
            '}',
            'td {',
            '    padding: 2px;',
            '    text-align: right;',
            '    border: 1px solid black;',
            '}',
            '</style>' ]
    else:
        headstart = ''
        headsep = ' '
        headend = ''
        colstart = ''
        colsep = ' '
        colend = ''
        css = list()
    #
    # Alias should be a dictionary that maps structure names to column names
    #
    if alias is None:
        #
        # Create a dummy alias dictionary
        #
        alias = dict(list(zip(array.dtype.names,array.dtype.names)))
    else:
        #
        # Fill in any missing values of the alias dictionary
        #
        for tag in array.dtype.names:
            if tag not in alias:
                alias[tag] = tag
    #
    # Formatcodes allows an override for certain columns.
    #
    if formatcodes is None:
        formatcodes = dict()
    #
    # This dictionary will hold the number of characters in each column
    #
    nchar = dict()
    #
    # Construct format codes for each column
    #
    for k, tag in enumerate(array.dtype.names):
        if tag in formatcodes:
            thiscode = formatcodes[tag]
            thisn = len(thiscode.format(array[tag][0]))
        else:
            d = array.dtype.fields[tag][0]
            if d.kind == 'i':
                thisn = max(max(len(str(array[tag].min())),len(str(array[tag].max()))),len(tag))
                thiscode = "{{{0:d}:{1:d}d}}".format(k,thisn)
            elif d.kind == 'f':
                if d.itemsize == 8:
                    prec = ddigit
                else:
                    prec = fdigit
                thisn = prec + 6
                if array[tag].min() < 0:
                    thisn += 1
                thiscode = "{{{0:d}:{1:d}.{2:d}g}}".format(k,thisn,prec)
            elif d.kind == 'S':
                thisn = max(d.itemsize,len(tag))
                thiscode = "{{{0:d}:{1:d}s}}".format(k,thisn)
            else:
                raise PydlutilsException("Unsupported kind: {0}".format(d.kind))
            formatcodes[tag] = thiscode
        nchar[tag] = thisn
    #
    # Start building an array of lines
    #
    lines = list()
    #
    # Construct header lines
    #
    if html:
        lines.append('<table>')
        hdr1 = headstart + headsep.join([alias[tag] for tag in array.dtype.names]) + headend
        lines.append(hdr1)
    else:
        if not no_head:
            hdr1 = headstart + headsep.join([("{{0:{0:d}s}}".format(nchar[tag])).format(alias[tag]) for tag in array.dtype.names]) + headend
            hdr2 = headstart + headsep.join(['-' * nchar[tag] for tag in array.dtype.names]) + headend
            lines.append(hdr1)
            lines.append(hdr2)
    #
    # Create a format string for the data from the individual format codes
    #
    rowformat = colstart + colsep.join([formatcodes[tag] for tag in array.dtype.names]) + colend
    if debug:
        print(rowformat)
    for k in range(array.size):
        lines.append(rowformat.format(*([decode_mixed(l) for l in array[k].tolist()])))
    if html:
        lines.append('</table>')
    if f is None:
        if not silent:
            print("\n".join(lines)+"\n")
    else:
        f.write("\n".join(lines)+"\n")
        f.close()
    return (lines,css)
