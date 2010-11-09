#
#
#
def djs_readcol(name,**kwargs):
    """Read a free-format ASCII file into pylab arrays.

    Lines that do not match the specified format are ignored
    (such as comments). Columns may be separated by commas or whitespace.
    Format codes:
    A = string
    B = byte
    D = double precision (float64)
    F = single precision (float32)
    I = integer (int16) (for compatibility with IDL)
    L = long (int32) (for compatibility with IDL)
    K = long64 (int64) (new!)
    X = skip
    """
    import re
    #
    # Number of lines
    #
    try:
        f = open(name,'r')
    except IOError:
        return None
    lines = f.readlines()
    f.close()
    nlines = len(lines)
    if 'silent' in kwargs:
        silent = True
    else:
        silent = False
    if 'debug' in kwargs:
        debug = True
    else:
        debug = False
    if debug:
        print "%s contains %d lines." % (name,nlines)
    if 'skip' in kwargs:
        skip = kwargs['skip']
    else:
        skip = 0
    nlines -= skip
    if 'numline' in kwargs:
        nlines = min(kwargs['numline'],nlines)
    #
    # Get the number of columns from the first non-skipped line
    #
    k = skip
    while lines[k][0] == '#':
        k += 1
    whitespace = re.compile(r'\s+')
    baseline = lines[k].strip().replace(',',' ')
    basecols = whitespace.split(baseline)
    ncol = len(basecols)
    if 'format' in kwargs:
        if re.match(r'^\(?[ABDFILX, ]+\)?$',kwargs['format'],re.IGNORECASE) is None:
            print "Invalid format string!"
            return None
        format = kwargs['format'].replace(' ','').upper().lstrip('(').rstrip(')').split(',')
        saveformat = [f for f in format if f != 'X']
        if len(format) < ncol:
            if not silent:
                print 'Format string has fewer columns than the file.'
            ncol = len(format)
    else:
        #
        # Assume all floating point format
        #
        format = list('F'*ncol)
        saveformat = format
    if debug:
        print ','.join(format)
    nread = 0
    goodlist = list()
    for l in lines[skip:nlines]:
        nread += 1
        if debug:
            print l
        if len(l) < ncol or l[0] == '#':
            if not silent:
                print 'Skipping line %d' % (skip+nread+1,)
            continue
        #
        # Split the line
        #
        cols = whitespace.split(l.strip().replace(',',' '))
        savecols = [cols[k] for k in range(ncol) if format[k] != 'X']
        savelist = list()
        if len(savecols) == len(saveformat):
            for k in range(len(saveformat)):
                if saveformat[k] == 'A':
                    #
                    # Save strings as is.
                    #
                    saved = savecols[k]
                elif saveformat[k] == 'B' or saveformat[k] == 'I' or saveformat[k] == 'L':
                    try:
                        saved = int(savecols[k])
                    except ValueError:
                        #
                        # Error, bad format, skip this line
                        #
                        if not silent:
                            print 'Skipping line %d' % skip+nread+1
                        continue
                elif saveformat[k] == 'F' or saveformat[k] == 'D':
                    try:
                        saved = float(savecols[k])
                    except ValueError:
                        #
                        # Error, bad format, skip this line
                        #
                        if not silent:
                            print 'Skipping line %d' % skip+nread+1
                        continue
                else:
                    print "Whoops, bad format! How did that happen?"
                    continue
                savelist.append(saved)
            if len(savelist) != len(saveformat):
                if not silent:
                    print "Skipping line %d" % skip+nread+1
        else:
            #
            # Error, not enough columns
            #
            if not silent:
                print "Skipping line %d" % skip+nread+1
            continue
        goodlist.append(savelist)
    if len(goodlist) == 0:
        raise IOError('No valid lines found for specified format')
    if not silent:
        print "%d valid lines read." % len(goodlist)
    #
    # Zip the good list
    #
    goodcols = zip(*goodlist)
    #
    # Convert the columns to pylab arrays
    #
    dtypes = { 'A':'S','B':'b','I':'i2','L':'i4','K':'i8','F':'f','D':'d' }
    converted = [np.array(goodcols[k],dtype=dtypes[saveformat[k]])
        for k in range(len(saveformat))]
    return tuple(converted)

