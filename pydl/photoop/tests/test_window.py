# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import os
import numpy as np
from astropy.tests.helper import raises
from .. import PhotoopException
from ..window import sdss_score, window_read, window_score


class TestWindow(object):
    """Test the functions in pydl.photoop.window.
    """

    def setup(self):
        pass

    def teardown(self):
        pass

    def test_sdss_score(self):
        pass

    def test_window_read_no_photo_resolve(self, monkeypatch):
        monkeypatch.delenv('PHOTO_RESOLVE', raising=False)
        with raises(PhotoopException) as e:
            window_read()
        assert (str(e.value) ==
                'You have not set the environment variable PHOTO_RESOLVE!')

    def test_window_read_all(self, monkeypatch, mocker):
        monkeypatch.setenv('PHOTO_RESOLVE', '/another/fake/directory')
        flist_size = 10
        hhh = mocker.MagicMock()
        hhh.data = {'SCORE': np.zeros((flist_size,), dtype=np.int16) }
        hh = mocker.MagicMock()
        hh[0] = None
        hh[1] = hhh
        h = mocker.patch('astropy.io.fits.open')
        h.return_value = hh
        r = window_read(flist=True, rescore=False, blist=True,
                        bcaps=True, balkans=True, findx=True, bindx=True)
        assert h.call_count == 5
        h.assert_any_call('/another/fake/directory/window_flist.fits')
        h.assert_any_call('/another/fake/directory/window_blist.fits')
        h.assert_any_call('/another/fake/directory/window_bcaps.fits')
        h.assert_any_call('/another/fake/directory/window_findx.fits')
        h.assert_any_call('/another/fake/directory/window_bindx.fits')

    def test_window_score_no_photo_calib(self, monkeypatch):
        monkeypatch.delenv('PHOTO_CALIB', raising=False)
        with raises(PhotoopException) as e:
            window_score()
        assert (str(e.value) ==
                'You have not set the environment variable PHOTO_CALIB!')

    def test_window_score_no_photo_resolve(self, monkeypatch):
        monkeypatch.setenv('PHOTO_CALIB', '/fake/directory')
        monkeypatch.delenv('PHOTO_RESOLVE', raising=False)
        with raises(PhotoopException) as e:
            window_score()
        assert (str(e.value) ==
                'You have not set the environment variable PHOTO_RESOLVE!')

    def test_window_score_no_window_flist(self, monkeypatch):
        monkeypatch.setenv('PHOTO_CALIB', '/fake/directory')
        monkeypatch.setenv('PHOTO_RESOLVE', '/another/fake/directory')
        with raises(PhotoopException) as e:
            window_score()
        assert str(e.value) == 'Unable to read FLIST file.'

    def test_window_score_fits_open(self, monkeypatch, mocker):
        monkeypatch.setenv('PHOTO_CALIB', '/fake/directory')
        monkeypatch.setenv('PHOTO_RESOLVE', '/another/fake/directory')
        score_size = 10
        hhh = mocker.MagicMock()
        hhh.data = {'SCORE': np.zeros((score_size,), dtype=np.int16) }
        hh = mocker.MagicMock()
        hh[0] = None
        hh[1] = hhh
        h = mocker.patch('astropy.io.fits.open')
        h.return_value = hh
        s = mocker.patch('pydl.photoop.window.sdss_score')
        s.return_value = np.zeros((score_size,), dtype=np.int16)
        window_score()
        h.assert_called_once_with('/another/fake/directory/window_flist.fits',
                                  mode='update')
        s.assert_called_once_with(h.return_value)

    def test_window_score_fits_open_rescore(self, monkeypatch, mocker):
        monkeypatch.setenv('PHOTO_CALIB', '/fake/directory')
        monkeypatch.setenv('PHOTO_RESOLVE', '/another/fake/directory')
        score_size = 10
        hhh = mocker.MagicMock()
        hhh.data = {'SCORE': np.zeros((score_size,), dtype=np.int16) }
        hh = mocker.MagicMock()
        hh[0] = None
        hh[1] = hhh
        h = mocker.patch('astropy.io.fits.open')
        h.return_value = hh
        s = mocker.patch('pydl.photoop.window.sdss_score')
        s.return_value = np.zeros((score_size,), dtype=np.int16)
        window_score(rescore=True)
        assert s.call_count == 1
        h.assert_called_once_with('/another/fake/directory/window_flist.fits',
                                  mode='readonly')
        s.assert_called_once_with(h.return_value)
