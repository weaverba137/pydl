# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test the functions in pydl.photoop.window.
"""
import os
import pytest
import numpy as np
from astropy.tests.helper import raises
from .. import PhotoopException
from ..window import sdss_score, window_read, window_score
from ...pydlutils.mangle import FITS_polygon


@pytest.fixture
def fits_open(request):
    data_size = 10
    m = request.getfixturevalue("mocker")
    hhh = m.MagicMock()
    hhh.data = {'SCORE': np.zeros((data_size,), dtype=np.int16),
                'ICAP': np.arange(data_size, dtype=np.int16),
                'NCAPS': np.ones((data_size,), dtype=np.int16),
                'X': np.zeros((data_size,), dtype=np.float32),
                'CM': np.zeros((data_size,), dtype=np.float32)}
    hh = m.MagicMock()
    hh.__enter__.return_value = [None, hhh]
    h = m.patch('astropy.io.fits.open')
    h.return_value = hh
    h.data_size = data_size
    return h


@pytest.fixture
def table_read(request):
    data_size = 10
    m = request.getfixturevalue("mocker")
    class fake_table(dict):
        def __len__(self):
            return data_size
    hh = fake_table()
    hh['SCORE'] = np.zeros((data_size,), dtype=np.int16)
    hh['IPRIMARY'] = np.zeros((data_size,), dtype=np.int16)
    hh['IBINDX'] = np.zeros((data_size,), dtype=np.int16)
    hh['ICAP'] = 4*np.arange(data_size, dtype=np.int16)
    hh['NCAPS'] = np.ones((data_size,), dtype=np.int16) + 3
    hh['WEIGHT'] = np.ones((data_size,), dtype=np.float32)
    hh['STR'] = np.zeros((data_size,), dtype=np.float32)
    hh['X'] = np.zeros((4*data_size, 3), dtype=np.float32)
    hh['CM'] = np.zeros((4*data_size,), dtype=np.float32)
    h = m.patch('astropy.table.Table.read')
    h.return_value = hh
    h.data_size = data_size
    return h


def test_sdss_score():
    pass


def test_window_read_no_photo_resolve(monkeypatch):
    monkeypatch.delenv('PHOTO_RESOLVE', raising=False)
    with raises(PhotoopException) as e:
        window_read()
    assert (str(e.value) ==
            'You have not set the environment variable PHOTO_RESOLVE!')


def test_window_read_all(monkeypatch, mocker, table_read):
    monkeypatch.setenv('PHOTO_RESOLVE', '/another/fake/directory')
    r = window_read(flist=True, rescore=False, blist=True,
                    bcaps=True, balkans=True, findx=True, bindx=True)
    assert table_read.call_count == 5
    table_read.assert_any_call('/another/fake/directory/window_flist.fits', hdu=1)
    table_read.assert_any_call('/another/fake/directory/window_blist.fits', hdu=1)
    table_read.assert_any_call('/another/fake/directory/window_bcaps.fits', hdu=1)
    table_read.assert_any_call('/another/fake/directory/window_findx.fits', hdu=1)
    table_read.assert_any_call('/another/fake/directory/window_bindx.fits', hdu=1)
    for k in ('flist', 'blist', 'bcaps', 'findx', 'bindx', 'balkans'):
        assert k in r
    for k in ('x', 'cm', 'use_caps'):
        assert hasattr(r['balkans'], k)


def test_window_read_balkans_only(monkeypatch, mocker, table_read):
    monkeypatch.setenv('PHOTO_RESOLVE', '/another/fake/directory')
    r = window_read(flist=True, rescore=False, balkans=True)
    assert table_read.call_count == 3
    table_read.assert_any_call('/another/fake/directory/window_flist.fits', hdu=1)
    table_read.assert_any_call('/another/fake/directory/window_blist.fits', hdu=1)
    table_read.assert_any_call('/another/fake/directory/window_bcaps.fits', hdu=1)
    for k in ('flist', 'balkans'):
        assert k in r
    assert isinstance(r['balkans'], FITS_polygon)


def test_window_read_rescore(monkeypatch, mocker, table_read):
    monkeypatch.setenv('PHOTO_RESOLVE', '/another/fake/directory')
    e = mocker.patch('os.path.exists')
    e.return_value = False
    w = mocker.patch('pydl.photoop.window.window_score')
    r = window_read(flist=True, rescore=True)
    assert table_read.call_count == 1
    table_read.assert_any_call('/another/fake/directory/window_flist_rescore.fits', hdu=1)
    e.assert_called_once_with('/another/fake/directory/window_flist_rescore.fits')
    w.assert_called_once_with(rescore=True)
    for k in ('flist',):
        assert k in r


def test_window_score_no_photo_calib(monkeypatch):
    monkeypatch.delenv('PHOTO_CALIB', raising=False)
    with raises(PhotoopException) as e:
        window_score()
    assert (str(e.value) ==
            'You have not set the environment variable PHOTO_CALIB!')


def test_window_score_no_photo_resolve(monkeypatch):
    monkeypatch.setenv('PHOTO_CALIB', '/fake/directory')
    monkeypatch.delenv('PHOTO_RESOLVE', raising=False)
    with raises(PhotoopException) as e:
        window_score()
    assert (str(e.value) ==
            'You have not set the environment variable PHOTO_RESOLVE!')


def test_window_score_no_window_flist(monkeypatch):
    monkeypatch.setenv('PHOTO_CALIB', '/fake/directory')
    monkeypatch.setenv('PHOTO_RESOLVE', '/another/fake/directory')
    with raises(PhotoopException) as e:
        window_score()
    assert str(e.value) == 'Unable to read FLIST file.'


def test_window_score_fits_open(monkeypatch, mocker, fits_open):
    monkeypatch.setenv('PHOTO_CALIB', '/fake/directory')
    monkeypatch.setenv('PHOTO_RESOLVE', '/another/fake/directory')
    s = mocker.patch('pydl.photoop.window.sdss_score')
    s.return_value = np.zeros((fits_open.data_size,), dtype=np.int16)
    window_score()
    fits_open.assert_called_once_with('/another/fake/directory/window_flist.fits',
                                      mode='update')
    s.assert_called_once_with(fits_open.return_value)


def test_window_score_fits_open_rescore(monkeypatch, mocker, fits_open):
    monkeypatch.setenv('PHOTO_CALIB', '/fake/directory')
    monkeypatch.setenv('PHOTO_RESOLVE', '/another/fake/directory')
    s = mocker.patch('pydl.photoop.window.sdss_score')
    s.return_value = np.zeros((fits_open.data_size,), dtype=np.int16)
    window_score(rescore=True)
    assert s.call_count == 1
    fits_open.assert_called_once_with('/another/fake/directory/window_flist.fits',
                                      mode='readonly')
    fits_open.return_value.writeto.assert_called_once_with('/another/fake/directory/window_flist_rescore.fits')
    s.assert_called_once_with(fits_open.return_value)
