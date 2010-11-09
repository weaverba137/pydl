#
#
#
def set_maskbits():
    """Populate the maskbits cache.
    """
    import yanny
    if pydlutils.sdss.maskbits is None:
        maskfile = yanny.yanny(os.getenv('IDLUTILS_DIR')+'/data/sdss/sdssMaskbits.par')
        if str(maskfile) == '':
            raise IOError('File with mask bits not found!')
        #
        # Parse the file & cache the results in maskbits
        #
        pydlutils.sdss.maskbits = dict()
        for k in range(maskfile.size('MASKBITS')):
            if maskfile['MASKBITS']['flag'][k] in pydlutils.sdss.maskbits:
                pydlutils.sdss.maskbits[maskfile['MASKBITS']['flag'][k]][maskfile['MASKBITS']['label'][k]] = maskfile['MASKBITS']['bit'][k]
            else:
                pydlutils.sdss.maskbits[maskfile['MASKBITS']['flag'][k]] = {
                maskfile['MASKBITS']['label'][k]:maskfile['MASKBITS']['bit'][k] }
    return


