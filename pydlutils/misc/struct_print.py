#
# $Id$
#
def struct_print(array,**kwargs):
    """Print a NumPy record array (analogous to an IDL structure) in a nice
    way.

    filename - a file name or file object
    """
    import numpy as np
    from pydlutils import PydlutilsException
    f = None # This variable will store a file handle
    if 'filename' in kwargs:
        if isinstance(kwargs['filename'],file):
            f = kwargs['filename']
        else:
            f = open(kwargs['filename'],'w')
    if 'fdigit' in kwargs:
        fdigit = kwargs['fdigit']
    else:
        fdigit = 5
    if 'ddigit' in kwargs:
        ddigit = kwargs['ddigit']
    else:
        ddigit = 7
    if 'html' in kwargs:
        html = True
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
        html = False
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
    if 'alias' in kwargs:
        if isinstance(kwargs['alias'],dict):
            alias = kwargs['alias']
            #
            # Fill in missing values of the alias dictionary
            #
            for tag in array.dtype.names:
                if tag not in alias:
                    alias[tag] = tag
        else:
            raise PydlutilsException("Invalid type for alias keyword.")
    else:
        #
        # Create a dummy alias dictionary
        #
        alias = dict(zip(array.dtype.names,array.dtype.names))
    #
    # Formatcodes allows an override for
    #
    if 'formatcodes' in kwargs:
        if isinstance(kwargs['formatcodes'],dict):
            formatcodes = kwargs['formatcodes']
        else:
            raise PydlutilsException("Invalid type for formatcodes keyword.")
    else:
        #
        # Create an empty dictionary that will hold the format codes.
        #
        formatcodes = dict()
    #
    # This dictionary will hold the number of characters in each column
    #
    nchar = dict()
    #
    # Construct format codes for each column
    #
    for tag in array.dtype.names:
        if tag in formatcodes:
            thiscode = formatcodes[tag]
            thisn = len(thiscode % array[tag][0])
        else:
            d = array.dtype.fields[tag][0]
            if d.kind == 'i':
                thisn = max(max(len(str(array[tag].min())),len(str(array[tag].max()))),len(tag))
                thiscode = "%%%dd" % thisn
            elif d.kind == 'f':
                if d.itemsize == 8:
                    prec = ddigit
                else:
                    prec = fdigit
                thisn = prec + 6
                if array[tag].min() < 0:
                    thisn += 1
                thiscode = "%%%d.%dg" % (thisn,prec)
            elif d.kind == 'S':
                thisn = max(d.itemsize,len(tag))
                thiscode = "%%%ds" % thisn
            else:
                raise PydlutilsException("Unsupported kind: %s" % d.kind)
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
        if 'no_head' not in kwargs:
            hdr1 = headstart + headsep.join([("%%%ds" % nchar[tag]) % alias[tag] for tag in array.dtype.names]) + headend
            hdr2 = headstart + headsep.join(['-' * nchar[tag] for tag in array.dtype.names]) + headend
            lines.append(hdr1)
            lines.append(hdr2)
    #
    # Create a format string for the data from the individual format codes
    #
    rowformat = colstart + colsep.join([formatcodes[tag] for tag in array.dtype.names]) + colend
    if 'debug' in kwargs:
        print rowformat
    lines += [rowformat % array[k].tolist() for k in range(array.size)]
    if html:
        lines.append('</table>')
    if f is None:
        print "\n".join(lines)
        print
    else:
        f.writelines(["%s\n" % l for l in lines])
        f.close()
    return (lines,css)
