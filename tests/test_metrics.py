"""
Tests for the `kpal.metrics` module.
"""


from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future import standard_library
from future.builtins import zip

import math

import numpy as np

from kpal import metrics

import utils

with standard_library.hooks():
    from collections import Counter


class TestMetrics(object):
    def test_distribution(self):
        a = np.random.randint(0, 21, 100)
        counts = Counter(a)
        assert metrics.distribution(a) == sorted(counts.items())

    def test_vector_length_float(self):
        a = np.random.rand(100)
        np.testing.assert_almost_equal(metrics.vector_length(a),
                                       math.sqrt(sum(x * x for x in a)))

    def test_vector_length_int(self):
        a = np.random.randint(0, 101, 100)
        np.testing.assert_almost_equal(metrics.vector_length(a),
                                       math.sqrt(sum(x * x for x in a)))

    def test_get_scale(self):
        a = np.random.randint(0, 101, 100)
        b = np.random.randint(0, 101, 100)

        scale_a, scale_b = metrics.get_scale(a, b)

        if a.sum() < b.sum():
            assert scale_a == b.sum() / a.sum()
            assert scale_b == 1.0
        else:
            assert scale_a == 1.0
            assert scale_b == a.sum() / b.sum()

    def test_scale_down(self):
        a = 1.0
        b = 1.0 + np.random.random()
        np.testing.assert_almost_equal(metrics.scale_down(a, b),
                                       (a / b, 1.0))

        a = 1.0 + np.random.random()
        b = 1.0
        np.testing.assert_almost_equal(metrics.scale_down(a, b),
                                       (1.0, b / a))

    def test_positive(self):
        a = np.random.randint(0, 21, 100)
        b = np.random.randint(0, 21, 100)

        np.testing.assert_array_equal(metrics.positive(a, b),
                                      [i if j else 0 for i, j in zip(a, b)])

    def test_multiset(self):
        a = np.random.randint(1, 101, 100)
        b = np.random.randint(1, 101, 100)
        pairwise = metrics.pairwise['prod']

        values = [pairwise(i, j) for i, j in zip(a, b) if i or j]
        np.testing.assert_almost_equal(metrics.multiset(a, b, pairwise),
                                       sum(values) / (len(values) + 1))

    def test_multiset_zeros(self):
        a = np.random.randint(0, 21, 100)
        b = np.random.randint(0, 21, 100)
        pairwise = metrics.pairwise['prod']

        values = [pairwise(i, j) for i, j in zip(a, b) if i or j]
        np.testing.assert_almost_equal(metrics.multiset(a, b, pairwise),
                                       sum(values) / (len(values) + 1))

    def test_multiset_many_zeros(self):
        a = np.random.randint(0, 3, 100)
        b = np.random.randint(0, 3, 100)
        pairwise = metrics.pairwise['prod']

        values = [pairwise(i, j) for i, j in zip(a, b) if i or j]
        np.testing.assert_almost_equal(metrics.multiset(a, b, pairwise),
                                       sum(values) / (len(values) + 1))

    def test_euclidean(self):
        a = np.random.randint(1, 101, 100)
        b = np.random.randint(1, 101, 100)

        np.testing.assert_almost_equal(metrics.euclidean(a, b),
                                       math.sqrt(sum((i - j) ** 2
                                                     for i, j in zip(a, b))))

    def test_cosine_similarity(self):
        a = np.random.randint(1, 101, 100)
        b = np.random.randint(1, 101, 100)

        cs = sum([x[0] * x[1] for x in zip(a, b)]) / (math.sqrt(sum(x * x for x in a)) *
                                                      math.sqrt(sum(x * x for x in b)))
        np.testing.assert_almost_equal(metrics.cosine_similarity(a, b),
                                       cs)

    def test_pairwise_prod(self):
        i = np.random.randint(100)
        j = np.random.randint(100)

        prod = abs(i - j) / float((i + 1) * (j + 1))
        np.testing.assert_almost_equal(metrics.pairwise['prod'](i, j), prod)

    def test_pairwise_sum(self):
        i = np.random.randint(100)
        j = np.random.randint(100)

        sum = abs(i - j) / float(i + j + 1)
        np.testing.assert_almost_equal(metrics.pairwise['sum'](i, j), sum)

    def test_summary_average(self):
        a = np.random.rand(100)
        np.testing.assert_almost_equal(metrics.summary['average'](a),
                                       a.mean())
