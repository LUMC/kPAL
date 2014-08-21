"""
Tests for the `k_mer.metrics` module.
"""


import numpy as np

from k_mer import metrics


class TestMetrics():
    def test_median_float_even(self):
        a = np.random.rand(100)
        assert metrics.median(a) == np.median(a)

    def test_median_float_odd(self):
        a = np.random.rand(101)
        assert metrics.median(a) == np.median(a)

    def test_median_int_even(self):
        a = np.random.random_integers(0, 100, 100)
        assert metrics.median(a) == np.median(a)

    def test_median_int_odd(self):
        a = np.random.random_integers(0, 100, 101)
        assert metrics.median(a) == np.median(a)
