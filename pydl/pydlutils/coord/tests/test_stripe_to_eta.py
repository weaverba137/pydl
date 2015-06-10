# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_stripe_to_eta():
    from .. import stripe_to_eta
    eta = stripe_to_eta(82)
    assert eta == -32.5
