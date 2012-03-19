#
#
#
def file_lines(filename,compress=False):
    """Replacement for IDL file_lines function.
    """
    if compress:
        import gzip
        with gzip.open(filename) as f:
            return len(f.readlines())
    else:
        with open(filename) as f:
            return len(f.readlines())
