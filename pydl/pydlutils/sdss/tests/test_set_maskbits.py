# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_set_maskbits():
    from ..set_maskbits import set_maskbits
    from os.path import dirname, join
    maskbits = set_maskbits(maskbits_file=join(dirname(__file__),'t','testMaskbits.par'))
    assert (set(maskbits.keys()) ==
        set(['TARGET', 'BOSS_TARGET1', 'PRIMTARGET', 'ANCILLARY_TARGET1',
            'ZWARNING', 'TTARGET', 'SECTARGET', 'LEGACY_TARGET2',
            'LEGACY_TARGET1', 'SPECIAL_TARGET2', 'FLUXMATCH_STATUS']))
    assert (set(maskbits['TARGET'].keys()) ==
        set(['QSO_FIRST_SKIRT', 'QSO_CAP', 'GALAXY_RED', 'STAR_CARBON',
            'STAR_WHITE_DWARF', 'GALAXY_RED_II', 'GALAXY_BIG',
            'GALAXY_BRIGHT_CORE', 'SERENDIP_MANUAL', 'STAR_SUB_DWARF',
            'QSO_FIRST_CAP', 'QSO_SKIRT', 'STAR_PN', 'STAR_BHB', 'QSO_HIZ',
            'STAR_BROWN_DWARF', 'SERENDIP_FIRST', 'SOUTHERN_SURVEY',
            'STAR_RED_DWARF', 'STAR_CATY_VAR', 'QSO_REJECT', 'GALAXY',
            'SERENDIP_RED', 'SERENDIP_DISTANT', 'QSO_MAG_OUTLIER', 'ROSAT_A',
            'ROSAT_C', 'ROSAT_B', 'ROSAT_E', 'ROSAT_D', 'SERENDIP_BLUE']))
