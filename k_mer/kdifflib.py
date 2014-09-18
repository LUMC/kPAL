"""
*k*-mer profile difference library.
"""


from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import str

import numpy as np

from . import metrics


class kMerDiff(object):
    """
    Class of distance functions.
    """
    def __init__(self, do_balance=False, do_positive=False, do_smooth=False,
                 summary=metrics.summary['min'], threshold=0, do_scale=False,
                 down=False, distance_function=None,
                 pairwise=metrics.pairwise['diff-prod']):
        """
        Initialise the class.

        :arg do_balance: Balance the profiles.
        :type do_balance: bool
        :arg do_positive: Only use positive values.
        :type do_positive: bool
        :arg do_smooth: Use dynamic smoothing.
        :type do_smooth: bool
        :arg summary: Summary function for dynamic smoothing.
        :type summary: function
        :arg threshold: Threshold for the summary function.
        :type threshold: int
        :arg do_scale: Scale the profiles.
        :type do_scale: bool
        :arg down: Normalise the scaling factors between 0 and 1.
        :type down: bool
        :arg distance_function: Use a specific distance function.
        :type distance_function: function
        :arg pairwise: Pairwise distance function for the multiset distance.
        :type pairwise: int
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

    def _dynamic_smooth(self, profile_left, profile_right, start, length):
        """
        Smooth two profiles by collapsing sub-profiles that do not meet the
        requirements governed by the selected summary function and the
        threshold.

        :arg profile_left, profile_right: A *k*-mer profile.
        :type profile_left, profile_right: klib.Profile
        :arg start: Start of the sub-profile to smooth.
        :type start: int
        :arg length: Length of the sub-profile to smooth.
        :type length: int
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

        c_left = self._collapse(profile_left.counts, start, length)
        c_right = self._collapse(profile_right.counts, start, length)

        if min(self._function(c_left),
               self._function(c_right)) <= self._threshold:
            # Collapse the sub-profile.
            profile_left.counts[start] = c_left.sum()
            profile_right.counts[start] = c_right.sum()

            # Remove the k-mer counts used to collapse.
            profile_left.counts[start + 1:start + length] = 0
            profile_right.counts[start + 1:start + length] = 0

            return

        new_length = length // 4
        for i in range(4):
            self._dynamic_smooth(profile_left, profile_right,
                                 start + i * new_length, new_length)

    def dynamic_smooth(self, profile_left, profile_right):
        """
        Smooth two profiles by collapsing sub-profiles that do not meet the
        requirements governed by the selected summary function and the
        threshold.

        :arg profile_left, profile_right: A *k*-mer profile.
        :type profile_left, profile_right: klib.Profile
        """
        self._dynamic_smooth(profile_left, profile_right, 0,
                             profile_left.number)

    def distance(self, profile_left, profile_right):
        """
        Calculate the distance between two *k*-mer profiles.

        :arg profile_left, profile_right: A *k*-mer profile.
        :type profile_left, profile_right: klib.Profile

        :return: The distance between `profile_left` and `profile_right`.
        :rtype: float
        """
        # No pun intended.
        copy_left = profile_left.copy()
        copy_right = profile_right.copy()

        if self._do_balance:
            copy_left.balance()
            copy_right.balance()

        if self._do_positive:
            copy_left.counts = metrics.positive(copy_left.counts,
                                                copy_right.counts)
            copy_right.counts = metrics.positive(copy_right.counts,
                                                 copy_left.counts)

        if self._do_smooth:
            self.dynamic_smooth(copy_left, copy_right)
        if self._do_scale:
            scale_left, scale_right = metrics.get_scale(copy_left.counts,
                                                        copy_right.counts)

            if self._down:
                scale_left, scale_right = metrics.scale_down(scale_left,
                                                             scale_right)
            copy_left.counts = copy_left.counts * scale_left
            copy_right.counts = copy_right.counts * scale_right

        if not self._distance_function:
            return metrics.multiset(copy_left.counts, copy_right.counts,
                                    self._pairwise)
        return self._distance_function(copy_left.counts, copy_right.counts)


def distance_matrix(profiles, output, precision, k_diff):
    """
    Make a distance matrix any number of *k*-mer profiles.

    :arg profiles: List of *k*-mer profiles.
    :type profiles: list(klib.Profile)
    :arg output: Open writable file handle.
    :type output: file-like object
    :arg precision: Number of digits in the output.
    :type precision: int
    :arg k_diff: A kMerDiff object.
    :type k_diff: kMerDiff
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
                         k_diff.distance(profiles[i], profiles[j])))

        output.write('\n')
