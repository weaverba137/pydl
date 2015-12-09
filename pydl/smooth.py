# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-


def smooth(signal, owidth, edge_truncate=False):
    """Replicates the IDL ``SMOOTH()`` function.

    Parameters
    ----------
    signal : array-like
        The array to be smoothed.
    owidth : :class:`int` or array-like
        Width of the smoothing window.  Can be a scalar or an array with
        length equal to the number of dimensions of `signal`.
    edge_truncate : :class:`bool`, optional
        Set `edge_truncate` to ``True`` to apply smoothing to all points.
        Points near the edge are normally excluded from smoothing.

    Returns
    -------
    array-like
        A smoothed array with the same dimesions and type as `signal`.

    Notes
    -----

    References
    ----------
    http://www.exelisvis.com/docs/SMOOTH.html

    Examples
    --------
    """
    if owidth % 2 == 0:
        width = owidth + 1
    else:
        width = owidth
    if width < 3:
        return signal
    n = signal.size
    istart = int((width-1)/2)
    iend = n - int((width+1)/2)
    w2 = int(width/2)
    s = signal.copy()
    for i in range(n):
        if i < istart:
            if edge_truncate:
                s[i] = (signal[0:istart+i+1].sum() +
                        (istart-i)*signal[0])/float(width)
        elif i > iend:
            if edge_truncate:
                s[i] = (signal[i-istart:n].sum() +
                        (i-iend)*signal[n-1])/float(width)
        else:
            s[i] = signal[i-w2:i+w2+1].sum()/float(width)
    return s
