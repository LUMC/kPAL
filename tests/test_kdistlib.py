"""
Tests for the `k_mer.metrics` module.
"""


from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future import standard_library

from io import StringIO

import numpy as np

from k_mer import kdistlib, klib, metrics

import utils

with standard_library.hooks():
    from collections import Counter


class TestKDistLib(object):
    def test_collapse(self):
        a = np.random.random_integers(0, 20, 100)
        start = 30
        length = 40

        step = length / 4
        expected = [sum(a[start + x[0]:start + x[1]])
                    for x in [(x * step, (x + 1) * step)
                              for x in range(4)]]

        k_dist = kdistlib.ProfileDistance()

        np.testing.assert_array_equal(k_dist._collapse(a, start, length),
                                      expected)

    def test_distance_matrix_one(self):
        counts = utils.counts(utils.SEQUENCES, 8)

        profiles = [klib.Profile(utils.as_array(counts, 8), 'a')]

        k_dist = kdistlib.ProfileDistance()
        out = StringIO()
        kdistlib.distance_matrix(profiles, out, 2, k_dist)

        assert out.getvalue().strip().split('\n') == ['1', 'a']

    def test_distance_matrix_two(self):
        counts_left = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_right = utils.counts(utils.SEQUENCES_RIGHT, 8)

        profiles = [klib.Profile(utils.as_array(counts_left, 8), 'a'),
                    klib.Profile(utils.as_array(counts_right, 8), 'b')]

        k_dist = kdistlib.ProfileDistance()
        out = StringIO()
        kdistlib.distance_matrix(profiles, out, 2, k_dist)

        assert out.getvalue().strip().split('\n') == ['2', 'a', 'b', '0.46']

    def test_distance_matrix_three(self):
        counts_left = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_right = utils.counts(utils.SEQUENCES_RIGHT, 8)

        profiles = [klib.Profile(utils.as_array(counts_left, 8), 'a'),
                    klib.Profile(utils.as_array(counts_right, 8), 'b'),
                    klib.Profile(utils.as_array(counts_left, 8), 'c')]

        k_dist = kdistlib.ProfileDistance()
        out = StringIO()
        kdistlib.distance_matrix(profiles, out, 2, k_dist)

        assert out.getvalue().strip().split('\n') == ['3', 'a', 'b', 'c', '0.46', '0.00 0.46']

    def test_ProfileDistance_dynamic_smooth(self):
        # If we use function=min and threshold=0, we should get the following
        # transformation:
        #
        #           | before           | after
        # ----------+------------------+-----------------
        #           | 0111111111111011 | 3000111111113000
        # profile A | ACGTACGTACGTACGT | ACGTACGTACGTACGT
        #           | AAAACCCCGGGGTTTT | AAAACCCCGGGGTTTT
        # ----------+------------------+-----------------
        #           | 0101111111111111 | 2000111111114000
        # profile B | ACGTACGTACGTACGT | ACGTACGTACGTACGT
        #           | AAAACCCCGGGGTTTT | AAAACCCCGGGGTTTT
        counts_a = Counter(['AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TG', 'TT'])
        counts_b = Counter(['AC', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT'])

        profile_a = klib.Profile(utils.as_array(counts_a, 2))
        profile_b = klib.Profile(utils.as_array(counts_b, 2))

        k_dist = kdistlib.ProfileDistance()
        k_dist.dynamic_smooth(profile_a, profile_b)

        counts_a = Counter(['AA', 'AA', 'AA', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TA', 'TA'])
        counts_b = Counter(['AA', 'AA', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TA', 'TA', 'TA'])

        np.testing.assert_array_equal(profile_a.counts, utils.as_array(counts_a, 2))
        np.testing.assert_array_equal(profile_b.counts, utils.as_array(counts_b, 2))

    def test_ProfileDistance_distance(self):
        counts_a = Counter(['AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TG', 'TT'])
        counts_b = Counter(['AC', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT'])

        profile_a = klib.Profile(utils.as_array(counts_a, 2))
        profile_b = klib.Profile(utils.as_array(counts_b, 2))

        k_dist = kdistlib.ProfileDistance()
        assert k_dist.distance(profile_a, profile_b) == 0.0625

    def test_ProfileDistance_distance_k8(self):
        counts_a = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_b = utils.counts(utils.SEQUENCES_RIGHT, 8)

        profile_a = klib.Profile(utils.as_array(counts_a, 8))
        profile_b = klib.Profile(utils.as_array(counts_b, 8))

        k_dist = kdistlib.ProfileDistance()
        np.testing.assert_almost_equal(k_dist.distance(profile_a, profile_b), 0.4626209322)

    def test_ProfileDistance_distance_unmodified(self):
        counts_a = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_b = utils.counts(utils.SEQUENCES_RIGHT, 8)

        profile_a = klib.Profile(utils.as_array(counts_a, 8))
        profile_b = klib.Profile(utils.as_array(counts_b, 8))

        k_dist = kdistlib.ProfileDistance(do_balance=True)
        k_dist.distance(profile_a, profile_b)

        utils.test_profile(profile_a, counts_a, 8)
        utils.test_profile(profile_b, counts_b, 8)
