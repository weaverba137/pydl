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
