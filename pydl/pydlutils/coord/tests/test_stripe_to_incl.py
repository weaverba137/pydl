# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_stripe_to_incl():
    from .. import stripe_to_incl
    incl = stripe_to_incl(82)
    assert incl == 0.0
    incl = stripe_to_incl(10)
    assert incl == 0.0
