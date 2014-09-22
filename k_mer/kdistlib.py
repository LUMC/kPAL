"""
*k*-mer profile distance library.

.. moduleauthor:: Leiden University Medical Center <humgen@lumc.nl>
.. moduleauthor:: Jeroen F.J. Laros <j.f.j.laros@lumc.nl>
.. moduleauthor:: Martijn Vermaat <m.vermaat.hg@lumc.nl>

.. Licensed under the MIT license, see the LICENSE file.
"""


from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import str

import numpy as np

from . import metrics


class ProfileDistance(object):
    """
    Class of distance functions.
    """
    def __init__(self, do_balance=False, do_positive=False, do_smooth=False,
                 summary=metrics.summary['min'], threshold=0, do_scale=False,
                 down=False, distance_function=None,
                 pairwise=metrics.pairwise['prod']):
        """
        Initialise the class.

        :arg bool do_balance: Balance the profiles.
        :arg bool do_positive: Only use positive values.
        :arg bool do_smooth: Use dynamic smoothing.
        :arg function summary: Summary function for dynamic smoothing.
        :arg int threshold: Threshold for the summary function.
        :arg bool do_scale: Scale the profiles.
        :arg bool down: Normalise the scaling factors between 0 and 1.
        :arg function distance_function: Use a specific distance function.
        :arg function pairwise: Pairwise distance function for the multiset
          distance.
        """
        self._do_balance = do_balance
        self._do_positive = do_positive
        self._do_smooth = do_smooth
        self._threshold = threshold
        self._do_scale = do_scale
        self._down = down
        self._distance_function = distance_function
        self._pairwise = pairwise
        self._function = summary

    def _collapse(self, vector, start, length):
        """
        Collapse a part of *k*-mer counts into a list of four numbers,
        representing the *k*-mers that start with a particular letter.

        :arg vector: Counts of a *k*-mer profile.
        :type vector: array_like, 1 dimension
        :arg start: Start of the area to collapse.
        :type start: int
        :arg length: Length of the area to collapse.
        :type length: int

        :return: Collapsed sub-profile.
        :rtype: ndarray
        """
        area = vector[start:start + length]
        return np.reshape(area, (4, length // 4)).sum(axis=1)

    def _dynamic_smooth(self, left, right, start, length):
        """
        Smooth two profiles by collapsing sub-profiles that do not meet the
        requirements governed by the selected summary function and the
        threshold.

        :arg k_mer.klib.Profile left, right: Profile to smooth.
        :arg int start: Start of the sub-profile to smooth.
        :arg int length: Length of the sub-profile to smooth.
        """
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
        if length == 1:
            return

        left_c = self._collapse(left.counts, start, length)
        right_c = self._collapse(right.counts, start, length)

        if min(self._function(left_c),
               self._function(right_c)) <= self._threshold:
            # Collapse the sub-profile.
            left.counts[start] = left_c.sum()
            right.counts[start] = right_c.sum()

            # Remove the k-mer counts used to collapse.
            left.counts[start + 1:start + length] = 0
            right.counts[start + 1:start + length] = 0

            return

        new_length = length // 4
        for i in range(4):
            self._dynamic_smooth(left, right, start + i * new_length,
                                 new_length)

    def dynamic_smooth(self, left, right):
        """
        Smooth two profiles by collapsing sub-profiles that do not meet the
        requirements governed by the selected summary function and the
        threshold.

        :arg k_mer.klib.Profile left, right: Profiles to smooth.
        """
        self._dynamic_smooth(left, right, 0, left.number)

    def distance(self, left, right):
        """
        Calculate the distance between two *k*-mer profiles.

        :arg k_mer.klib.Profile left, right: Profiles to calculate distance
          between.

        :return: The distance between `left` and `right`.
        :rtype: float
        """
        left = left.copy()
        right = right.copy()

        if self._do_balance:
            left.balance()
            right.balance()

        if self._do_positive:
            left.counts = metrics.positive(left.counts, right.counts)
            right.counts = metrics.positive(right.counts, left.counts)

        if self._do_smooth:
            self.dynamic_smooth(left, right)
        if self._do_scale:
            left_scale, right_scale = metrics.get_scale(left.counts,
                                                        right.counts)

            if self._down:
                left_scale, right_scale = metrics.scale_down(left_scale,
                                                             right_scale)
            left.counts = left.counts * left_scale
            right.counts = right.counts * right_scale

        if not self._distance_function:
            return metrics.multiset(left.counts, right.counts, self._pairwise)
        return self._distance_function(left.counts, right.counts)


def distance_matrix(profiles, output, precision, dist):
    """
    Make a distance matrix for any number of *k*-mer profiles.

    :arg list(Profile) profiles: List of profiles.
    :arg output: Open writable file handle.
    :type output: file-like object
    :arg int precision: Number of digits in the output.
    :arg k_mer.kdistlib.ProfileDistance dist: A distance functions object.
    """
    input_count = len(profiles)

    print(str(input_count), file=output)
    for i in profiles:
        print(i.name, file=output)
    for i in range(1, input_count):
        for j in range(i):
            if (j):
                output.write(' ')
            output.write('{{0:.{0}f}}'.format(precision).format(
                         dist.distance(profiles[i], profiles[j])))

        output.write('\n')
