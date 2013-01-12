#!/usr/bin/python

"""
k-mer profile difference library.
"""

import kLib
import metrics
import copy

class kMerDiff():
    """
    Class of distance functions.
    """
    def __init__(self, balance=False, positive=False, smooth=False, summary=0,
        threshold=0, scale=False, scaleDown=False, multiset=True, pairwise=0):
        """
        Initialise the class.

        @arg balance: Balance the profiles.
        @type balance: bool
        @arg positive: Only use positive values.
        @type positive: bool
        @arg smooth: Use dynamic smoothing.
        @type smooth: bool
        @arg summary: Summary function for dynamic smoothing.
        @type summary: int
        @arg threshold: Threshold for the summary function.
        @type threshold: int
        @arg scale: Scale the profiles.
        @type scale: bool
        @arg scaleDown: Normalise the scaling factors between 0 and 1.
        @type scaleDown: bool
        @arg multiset: Use the multiset distance metric, euclidean otherwise.
        @type multiset: bool
        @arg pairwise: Pairwise distance function for the multiset distance.
        @type pairwise: int
        """

        self.__balance = balance
        self.__positive = positive
        self.__smooth = smooth
        self.__threshold = threshold
        self.__scale = scale
        self.__scaleDown = scaleDown
        self.__multiset = multiset

        if pairwise not in range(len(metrics.pairwise)):
            raise ValueError("Invalid pairwise distance algorithm.")
        if summary not in range(len(metrics.summary)):
            raise ValueError("Invalid summary function.")

        self.__pairwise = metrics.pairwise[pairwise]
        self.__function = metrics.summary[summary]
    #__init__

    def __collapse(self, vector, start, length):
        """
        Collapse a part of k-mer counts into a list of four numbers,
        representing the k-mers that start with a particular letter.

        @arg vector: Counts of a k-mer profile.
        @type vector: list[float]
        @arg start: Start of the area to collapse.
        @type start: int
        @arg length: Length of the area to collapse.
        @type length: int

        @returns:
        @rtype: list[float]
        """
        step = length / 4

        return map(lambda x: sum(vector[start + x[0]:start + x[1]]),
            map(lambda x: (x * step, (x + 1) * step), range(4)))
    #__collapse

    def __dynamicSmooth(self, profile1, profile2, start, length):
        """
        Smooth two profiles by collapsing sub-profiles that do not meet the
        requirements governed by the selected summary function and the
        threshold.

        @arg profile1: A k-mer profile.
        @type profile1: object(kMer)
        @arg profile2: A k-mer profile.
        @type profile2: object(kMer)
        @arg start: Start of the sub-profile to smooth.
        @type start: int
        @arg length: Length of the sub-profile to smooth.
        @type length: int
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
        newLength = length / 4
        for i in range(4):
            self.__dynamicSmooth(profile1, profile2, start + i * newLength,
                newLength)
    #__dynamicSmooth

    def dynamicSmooth(self, profile1, profile2):
        """
        Smooth two profiles by collapsing sub-profiles that do not meet the
        requirements governed by the selected summary function and the
        threshold.

        @arg profile1: A k-mer profile.
        @type profile1: object(kMer)
        @arg profile2: A k-mer profile.
        @type profile2: object(kMer)
        """
        self.__dynamicSmooth(profile1, profile2, 0, profile1.number)
    #dynamicSmooth

    def calcDistance(self, profile1, profile2):
        """
        Calculate the distance between two k-mer profiles.

        @arg profile1: A k-mer profile.
        @type profile1: object(kMer)
        @arg profile2: A k-mer profile.
        @type profile2: object(kMer)

        @returns: The distance between {profile1} and {profile2}.
        @rtype: float
        """
        temp1 = copy.deepcopy(profile1)
        temp2 = copy.deepcopy(profile2)

        if self.__balance:
            temp1.balance()
            temp2.balance()
        #if
        if self.__positive:
            temp1.count = metrics.positive(temp1.count, temp2.count)
            temp2.count = metrics.positive(temp2.count, temp1.count)
        #if
        if self.__smooth:
            self.dynamicSmooth(temp1, temp2)
        if self.__scale:
            scale1, scale2 = metrics.calcScale(temp1.count, temp2.count)

            if self.__scaleDown:
                scale1, scale2 = metrics.scaleDown(scale1, scale2)
            temp1.count = metrics.scale(temp1.count, scale1)
            temp2.count = metrics.scale(temp2.count, scale2)
        #if
        if not self.__multiset:
            return metrics.euclideanDistance(temp1.count, temp2.count)
        return metrics.multisetDistance(temp1.count, temp2.count,
            self.__pairwise)
    #calcDistance
#kMerDiff

def makeDistanceMatrix(profiles, output, precision, kDiff):
    """
    Make a distance matrix any number of k-mer profiles.

    @arg profiles: List of k-mer profiles.
    @type profiles: list[kMer]
    @arg output: Open handle to a writable file.
    @type output: stream
    @arg precision: Number of digits in the output.
    @type precision: int
    @arg kDiff: A kMerDiff object.
    @type kDiff: object(kMerDiff)
    """
    numberOfInputs = len(profiles)

    output.write("%i\n" % numberOfInputs)
    for i in profiles:
        output.write("%s\n" % i.name)
    for i in range(1, numberOfInputs):
        for j in range(i):
            if (j):
                output.write(' ')
            output.write(("%%.%if" % precision) % 
                kDiff.calcDistance(profiles[i], profiles[j]))
        #for
        output.write('\n')
    #for
#makeDistanceMatrix
