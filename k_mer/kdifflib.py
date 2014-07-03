#!/usr/bin/python

"""
k-mer profile difference library.
"""

import copy

from . import klib, metrics

class kMerDiff():
    """
    Class of distance functions.
    """
    def __init__(self, do_balance=False, do_positive=False, do_smooth=False,
            summary=metrics.summary["min"], threshold=0, do_scale=False,
            down=False, distance_function=None,
            pairwise=metrics.pairwise["diff-prod"]):
        """
        Initialise the class.

        :arg do_balance: Balance the profiles.
        :type do_balance: bool
        :arg do_positive: Only use positive values.
        :type do_positive: bool
        :arg do_smooth: Use dynamic smoothing.
        :type do_smooth: bool
        :arg summary: Summary function for dynamic smoothing.
        :type summary: int
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

        self.__do_balance = do_balance
        self.__do_positive = do_positive
        self.__do_smooth = do_smooth
        self.__threshold = threshold
        self.__do_scale = do_scale
        self.__down = down
        self.__distance_function = distance_function
        self.__pairwise = pairwise
        self.__function = summary
    #__init__

    def __collapse(self, vector, start, length):
        """
        Collapse a part of k-mer counts into a list of four numbers,
        representing the k-mers that start with a particular letter.

        :arg vector: Counts of a k-mer profile.
        :type vector: list[float]
        :arg start: Start of the area to collapse.
        :type start: int
        :arg length: Length of the area to collapse.
        :type length: int

        :return: Collapsed sub-profile.
        :rtype: list[float]
        """
        step = length / 4

        return map(lambda x: sum(vector[start + x[0]:start + x[1]]),
            map(lambda x: (x * step, (x + 1) * step), range(4)))
    #__collapse

    def __dynamic_smooth(self, profile1, profile2, start, length):
        """
        Smooth two profiles by collapsing sub-profiles that do not meet the
        requirements governed by the selected summary function and the
        threshold.

        :arg profile1: A k-mer profile.
        :type profile1: object(kMer)
        :arg profile2: A k-mer profile.
        :type profile2: object(kMer)
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

        c1 = self.__collapse(profile1.count, start, length)
        c2 = self.__collapse(profile2.count, start, length)

        if min(self.__function(c1), self.__function(c2)) <= self.__threshold:
            # Collapse the sub-profile.
            profile1.count[start] = sum(c1)
            profile2.count[start] = sum(c2)

            # Remove the k-mer counts used to collapse.
            for i in range(start + 1, start + length):
                profile1.count[i] = 0
                profile2.count[i] = 0
            #for
            return
        #if
        new_length = length / 4
        for i in range(4):
            self.__dynamic_smooth(profile1, profile2, start + i * new_length,
                new_length)
    #__dynamic_smooth

    def dynamic_smooth(self, profile1, profile2):
        """
        Smooth two profiles by collapsing sub-profiles that do not meet the
        requirements governed by the selected summary function and the
        threshold.

        :arg profile1: A k-mer profile.
        :type profile1: object(kMer)
        :arg profile2: A k-mer profile.
        :type profile2: object(kMer)
        """
        self.__dynamic_smooth(profile1, profile2, 0, profile1.number)
    #dynamic_smooth

    def distance(self, profile1, profile2):
        """
        Calculate the distance between two k-mer profiles.

        :arg profile1: A k-mer profile.
        :type profile1: object(kMer)
        :arg profile2: A k-mer profile.
        :type profile2: object(kMer)

        :return: The distance between {profile1} and {profile2}.
        :rtype: float
        """
        temp1 = copy.deepcopy(profile1)
        temp2 = copy.deepcopy(profile2)

        if self.__do_balance:
            temp1.balance()
            temp2.balance()
        #if
        if self.__do_positive:
            temp1.count = metrics.positive(temp1.count, temp2.count)
            temp2.count = metrics.positive(temp2.count, temp1.count)
        #if
        if self.__do_smooth:
            self.dynamic_smooth(temp1, temp2)
        if self.__do_scale:
            scale1, scale2 = metrics.get_scale(temp1.count, temp2.count)

            if self.__down:
                scale1, scale2 = metrics.scale_down(scale1, scale2)
            temp1.count = metrics.scale(temp1.count, scale1)
            temp2.count = metrics.scale(temp2.count, scale2)
        #if
        if not self.__distance_function:
            return metrics.multiset(temp1.count, temp2.count,
                self.__pairwise)
        return self.__distance_function(temp1.count, temp2.count)
    #distance
#kMerDiff

def distance_matrix(profiles, output, precision, k_diff):
    """
    Make a distance matrix any number of k-mer profiles.

    :arg profiles: List of k-mer profiles.
    :type profiles: list[kMer]
    :arg output: Open handle to a writable file.
    :type output: stream
    :arg precision: Number of digits in the output.
    :type precision: int
    :arg k_diff: A kMerDiff object.
    :type k_diff: object(kMerDiff)
    """
    input_count = len(profiles)

    output.write("%i\n" % input_count)
    for i in profiles:
        output.write("%s\n" % i.name)
    for i in range(1, input_count):
        for j in range(i):
            if (j):
                output.write(' ')
            output.write(("%%.%if" % precision) % 
                k_diff.distance(profiles[i], profiles[j]))
        #for
        output.write('\n')
    #for
#distance_matrix
