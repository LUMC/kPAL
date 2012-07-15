#!/usr/bin/python

"""
k-mer profile difference library.
"""

import kLib
import math

class kMerDiff():
    """
    Class of distance functions.
    """
    algorithmHelp = "Distance algorithm to use (0 = multiset, " + \
        "1 = euclidean, 2 = positive multiset, 3 = relative multiset)."
    downHelp = "Scale down."
    algorithmError = "Invalid algorithm."

    def __init__(self, algorithm, down=False):
        """
        Initialise the class.

        @arg algorithm: Select which distance algorithm to use.
        @type algorithm: integer
        """
        self.distance = None
        self.down = down
        algorithms = [
            self.__multisetDistance,
            self.__euclideanDistance,
            self.__positiveMultisetDistance,
            self.__relativeMultisetDistance
        ]

        if algorithm not in range(len(algorithms)):
            return None

        self.distance = algorithms[algorithm]
    #__init__

    def __scale(self, profile1, profile2):
        """
        Calculate scaling factors based upon total counts. One of the factors
        is always one (the other is either one or larger than one).

        @arg profile1: A k-mer profile.
        @type profile1: object(kMer)
        @arg profile2: A k-mer profile.
        @type profile2: object(kMer)

        @returns: A tuple of scaling factors.
        @rtype: tuple(float)
        """
        scale1 = 1.0
        scale2 = 1.0

        # Calculate the scaling factors in such a way that no element is
        # between 0 and 1.
        if profile1.total < profile2.total:
            scale1 = float(profile2.total) / profile1.total
        else:
            scale2 = float(profile1.total) / profile2.total

        if self.down:
            factor = max(scale1, scale2)
            return scale1 / factor, scale2 / factor
        #if

        return scale1, scale2
    #__scale

    def __positiveScale(self, profile1, profile2):
        """
        Calculate scaling factors based upon counts that are non-negative in
        both k-mer sets. One of the factors is always one (the other is either
        one or larger than one).

        @arg profile1: A k-mer profile.
        @type profile1: object(kMer)
        @arg profile2: A k-mer profile.
        @type profile2: object(kMer)

        @returns: A tuple of scaling factors.
        @rtype: tuple(float)
        """
        scale1 = 1.0
        scale2 = 1.0
        total1 = 0
        total2 = 0

        if profile1.nonZero == profile1.number and \
            profile2.nonZero == profile2.number:
            return self.__scale(profile1, profile2)

        for i in range(len(profile1.count)):
            if profile1.count[i] and profile2.count[i]:
                total1 += profile1.count[i]
                total2 += profile2.count[i]
            #if

        # Calculate the scaling factors in such a way that no element is
        # between 0 and 1.
        if profile1.total < profile2.total:
            scale1 = float(total2) / total1
        else:
            scale2 = float(total1) / total2

        if self.down:
            factor = max(scale1, scale2)
            return scale1 / factor, scale2 / factor
        #if

        return scale1, scale2
    #__positiveScale

    def __collapse(self, l, start, length):
        """
        """
        step = length / 4

        return map(lambda x: sum(l[start + x[0]:start + x[1]]),
            map(lambda x: (x * step, (x + 1) * step), range(4)))
    #__collapse

    def dynamicSmooth(self, profile1, profile2, start, length, function,
        threshold):
        """
        Start with:
        profile1, profile2, 0, 4 ** profile1.length, function, threshold
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

        if min(function(c1), function(c2)) <= threshold:
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
            self.dynamicSmooth(profile1, profile2, start + i * newLength,
                newLength, function, threshold)
    #dynamicSmooth

    def __multisetDistance(self, profile1, profile2):
        """
        Calculate the multiset distance between two k-mer profiles.

        @arg profile1: A k-mer profile.
        @type profile1: object(kMer)
        @arg profile2: A k-mer profile.
        @type profile2: object(kMer)

        @returns: The multiset distance between {profile1} and {profile2}.
        @rtype: float
        """
        c = 0.0
        d = 1

        scale1, scale2 = self.__scale(profile1, profile2)

        # Calculate the counter and the denominator of the distance function.
        for i in range(profile1.number):
            if profile1.count[i] or profile2.count[i]:
                c += (abs(round(scale1 * profile1.count[i]) -
                          round(scale2 * profile2.count[i])) /
                    ((profile1.count[i] + 1) * (profile2.count[i] + 1)))
                d += 1
            #if
        #for

        return c / d
    #__multisetDistance

    def __euclideanDistance(self, profile1, profile2):
        """
        Calculate the Euclidean distance between two k-mer profiles.

        @arg profile1: A k-mer profile.
        @type profile1: object(kMer)
        @arg profile2: A k-mer profile.
        @type profile2: object(kMer)

        @returns: The Euclidean distance between {profile1} and {profile2}.
        @rtype: float
        """
        sumOfSquares = 1

        scale1, scale2 = self.__scale(profile1, profile2)

        # Calculate the counter and the denominator of the distance function.
        for i in range(profile1.number):
            sumOfSquares += (round(scale1 * profile1.count[i]) -
                round(scale2 * profile2.count[i])) ** 2

        return math.sqrt(sumOfSquares)
    #__euclideanDistance

    def __positiveMultisetDistance(self, profile1, profile2):
        """
        Calculate the positive multiset distance between two k-mer profiles.

        @arg profile1: A k-mer profile.
        @type profile1: object(kMer)
        @arg profile2: A k-mer profile.
        @type profile2: object(kMer)

        @returns: The positive multiset distance between {profile1} and
            {profile2}.
        @rtype: float
        """
        c = 0.0
        d = 1

        scale1, scale2 = self.__positiveScale(profile1, profile2)

        # Calculate the counter and the denominator of the distance function.
        for i in range(profile1.number):
            if profile1.count[i] and profile2.count[i]:
                c += (abs(round(scale1 * profile1.count[i]) -
                          round(scale2 * profile2.count[i])) /
                    ((profile1.count[i] + 1) * (profile2.count[i] + 1)))
                d += 1
            #if
        #for

        return c / d
    #__positiveMultisetDistance

    def __relativeMultisetDistance(self, profile1, profile2):
        """
        Calculate the relative multiset distance between two k-mer profiles.

        @arg profile1: A k-mer profile.
        @type profile1: object(kMer)
        @arg profile2: A k-mer profile.
        @type profile2: object(kMer)

        @returns: The relative multiset distance between {profile1} and
            {profile2}.
        @rtype: float
        """
        c = 0.0
        d = 1

        scale1, scale2 = self.__scale(profile1, profile2)

        for i in range(profile1.number):
            for j in range(i):
                diff1 = abs(profile1.count[i] - profile1.count[j])
                diff2 = abs(profile2.count[i] - profile2.count[j])

                if diff1 or diff2:
                    c += (abs(round(scale1 * diff1) - round(scale2 * diff2)) / 
                        ((diff1 + 1) * (diff2 + 1)))
                    d += 1
                #if
            #for

        return c / d
    #__relativeMultisetDistance
#kMerDiff

def makeDistanceMatrix(counts, output, precision, kDiff):
    """
    Make a distance matrix any number of k-mer profiles.

    @arg counts: List of k-mer profiles.
    @type counts: list[kMer]
    @arg output: Open handle to a writable file.
    @type output: stream
    @arg precision: Number of digits in the output.
    @type precision: integer
    @arg kDiff: A kMerDiff object.
    @type kDiff: list[kMerDiff]
    """
    numberOfInputs = len(counts)

    output.write("%i\n" % numberOfInputs)
    for i in range(1, numberOfInputs + 1):
        output.write("%i\n" % i)
    for i in range(1, numberOfInputs):
        for j in range(i):
            if (j):
                output.write(' ')
            output.write(("%%.%if" % precision) % 
                kDiff.distance(counts[i], counts[j]))
        #for
        output.write('\n')
    #for
#makeDistanceMatrix
