#
#
#
def smooth(signal,owidth,edge_truncate=False):
    """Replacement for IDL smooth function.
    """
    if owidth % 2 == 0:
       width = owidth + 1
    else:
       width = owidth
    if width < 3:
        return signal
    n = signal.size
    istart = (width-1)/2
    iend = n - (width+1)/2
    w2 = width/2
    s = signal.copy()
    for i in range(n):
        if i < istart:
            if edge_truncate:
                s[i] = (signal[0:istart+i+1].sum() + (istart-i)*signal[0])/float(width)
        elif i > iend:
            if edge_truncate:
                s[i] = (signal[i-istart:n].sum() + (i-iend)*signal[n-1])/float(width)
        else:
            s[i] = signal[i-w2:i+w2+1].sum()/float(width)
    return s

