# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
This subpackage corresponds to the sdssio directory of photoop.
"""
from .filtername import filtername
from .filternum import filternum
from .sdss_calib import sdss_calib
from .sdss_name import sdss_name
from .sdss_path import sdss_path
from .sdssflux2ab import sdssflux2ab
#
# Filename formats used by sdss_name and sdss_path
#
_name_formats = {
    'apObj':"{ftype}-{run:06d}-{filter}{camcol:1d}-{field:04d}.fit",
    'calibMatch':"{ftype}-{run:06d}-{camcol:1d}.fits",
    'calibPhotom':"{ftype}-{run:06d}-{camcol:1d}.fits",
    'calibPhotomGlobal':"{ftype}-{run:06d}-{camcol:1d}.fits",
    'fakeIdR':"idR-{run:06d}-{filter}{camcol:1d}-{field:04d}.fit",
    'fpAtlas':"{ftype}-{run:06d}-{camcol:1d}-{field:04d}.fit",
    'fpBIN':"{ftype}-{run:06d}-{filter}{camcol:1d}-{field:04d}.fit",
    'fpC':"{ftype}-{run:06d}-{filter}{camcol:1d}-{field:04d}.fit",
    'fpFieldStat':"{ftype}-{run:06d}-{camcol:1d}-{field:04d}.fit",
    'fpM':"{ftype}-{run:06d}-{filter}{camcol:1d}-{field:04d}.fit",
    'fpObjc':"{ftype}-{run:06d}-{camcol:1d}-{field:04d}.fit",
    'hoggObj':"{ftype}-{run:06d}-{camcol:1d}-{field:04d}.fits",
    'idFF':"{ftype}-{run:06d}-{filter}{camcol:1d}.fit",
    'idR':"{ftype}-{run:06d}-{filter}{camcol:1d}-{field:04d}.fit",
    'idRR':"{ftype}-{run:06d}-{filter}{camcol:1d}-{field:04d}.fit",
    'psBB':"{ftype}-{run:06d}-{filter}{camcol:1d}-{field:04d}.fit",
    'psFF':"{ftype}-{run:06d}-{filter}{camcol:1d}.fit",
    'psField':"{ftype}-{run:06d}-{camcol:1d}-{field:04d}.fit",
    'reObjGlobal':"{ftype}-{run:06d}-{camcol:1d}-{field:04d}.fits",
    'reObjRun':"{ftype}-{run:06d}-{camcol:1d}-{field:04d}.fits",
    'reObjTmp':"{ftype}-{run:06d}-{camcol:1d}-{field:04d}.fits",
    'tsField':"{ftype}-{run:06d}-{camcol:1d}-{rerun}-{field:04d}.fit",
    }
#
#
#
_path_formats = {
    'apObj':"{redux}/{rerun}/{run:d}/objcs/{camcol:1d}",
    'calibMatch':"{redux}/{rerun}/{run:d}/nfcalib",
    'calibPhotom':"{redux}/{rerun}/{run:d}/nfcalib",
    'calibPhotomGlobal':"{calib}/{rerun}/{run:d}/nfcalib",
    'fakeIdR':"{data}/{run:d}/fake_fields/{camcol:1d}",
    'fpAtlas':"{redux}/{rerun}/{run:d}/objcs/{camcol:1d}",
    'fpBIN':"{redux}/{rerun}/{run:d}/objcs/{camcol:1d}",
    'fpC':"{redux}/{rerun}/{run:d}/objcs/{camcol:1d}",
    'fpFieldStat':"{redux}/{rerun}/{run:d}/objcs/{camcol:1d}",
    'fpM':"{redux}/{rerun}/{run:d}/objcs/{camcol:1d}",
    'fpObjc':"{redux}/{rerun}/{run:d}/objcs/{camcol:1d}",
    'hoggObj':"{redux}/{rerun}/{run:d}/objcs/{camcol:1d}",
    'idFF':"{redux}/{rerun}/{run:d}/objcs/{camcol:1d}",
    'idR':"{data}/{run:d}/fields/{camcol:1d}",
    'idRR':"{data}/{run:d}/fields/{camcol:1d}",
    'psBB':"{redux}/{rerun}/{run:d}/objcs/{camcol:1d}",
    'psFF':"{redux}/{rerun}/{run:d}/objcs/{camcol:1d}",
    'psField':"{redux}/{rerun}/{run:d}/objcs/{camcol:1d}",
    'reObjGlobal':"{resolve}/{rerun}/{run:d}/resolve/{camcol:1d}",
    'reObjRun':"{redux}/{rerun}/{run:d}/resolve/{camcol:1d}",
    'reObjTmp':"{resolve}/{rerun}/{run:d}/resolve/{camcol:1d}",
    'tsField':"{redux}/{rerun}/{run:d}/calibChunks/{camcol:1d}",
    }
