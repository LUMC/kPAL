#!/usr/bin/python

"""
@requires: argparse
@requires: kMer
@requires: math
"""

import argparse # ArgumentParser().
import kMer     # kMer().
import math     # sqrt().

class kMerDiff() :
    """
    Class of distance functions.
    """

    algorithmHelp = "Distance algorithm to use (0 = multiset, " + \
        "1 = euclidean, 2 = positive multiset, 3 = relative multiset)."
    algorithmError = "Invalid algorithm."

    def __init__(self, algorithm) :
        """
        Initialise the class.

        @arg algorithm: Select which distance algorithm to use.
        @type algorithm: integer
        """

        self.distance = None
        algorithms = [
            self.__multisetDistance,
            self.__euclideanDistance,
            self.__positiveMultisetDistance,
            self.__relativeMultisetDistance
        ]

        if algorithm not in range(len(algorithms)) :
            return None

        self.distance = algorithms[algorithm]
    #__init__

    def __scale(self, kMerIn1, kMerIn2) :
        """
        Calculate scaling factors based upon total counts. One of the factors
        is always one (the other is either one or larger than one).

        @arg kMerIn1: A kMer instance.
        @type kMerIn1: object(kMer)
        @arg kMerIn2: A kMer instance.
        @type kMerIn2: object(kMer)

        @returns: A tuple of scaling factors.
        @rtype: tuple(float)
        """

        scale1 = 1.0
        scale2 = 1.0

        # Calculate the scaling factors in such a way that no element is
        # between 0 and 1.
        if kMerIn1.totalKMers < kMerIn2.totalKMers :
            scale1 = float(kMerIn2.totalKMers) / kMerIn1.totalKMers
        else :
            scale2 = float(kMerIn1.totalKMers) / kMerIn2.totalKMers

        return scale1, scale2
    #__scale

    def __positiveScale(self, kMerIn1, kMerIn2) :
        """
        Calculate scaling factors based upon counts that are non-negative in
        both k-mer sets. One of the factors is always one (the other is either
        one or larger than one).

        @arg kMerIn1: A kMer instance.
        @type kMerIn1: object(kMer)
        @arg kMerIn2: A kMer instance.
        @type kMerIn2: object(kMer)

        @returns: A tuple of scaling factors.
        @rtype: tuple(float)
        """

        scale1 = 1.0
        scale2 = 1.0
        total1 = 0
        total2 = 0

        if kMerIn1.nonZeroKMers == kMerIn1.numberOfKMers and \
            kMerIn2.nonZeroKMers == kMerIn2.numberOfKMers :
            return self.__scale(kMerIn1, kMerIn2)

        for i in range(len(kMerIn1.kMerCount)) :
            if kMerIn1.kMerCount[i] and kMerIn2.kMerCount[i] :
                total1 += kMerIn1.kMerCount[i]
                total2 += kMerIn2.kMerCount[i]
            #if

        # Calculate the scaling factors in such a way that no element is
        # between 0 and 1.
        if kMerIn1.totalKMers < kMerIn2.totalKMers :
            scale1 = float(total2) / total1
        else :
            scale2 = float(total1) / total2

        return scale1, scale2
    #__positiveScale

    def __multisetDistance(self, kMerIn1, kMerIn2) :
        """
        Calculate the multiset distance between two kMer instances.

        @arg kMerIn1: A kMer instance.
        @type kMerIn1: object(kMer)
        @arg kMerIn2: A kMer instance.
        @type kMerIn2: object(kMer)

        @returns: The multiset distance between {kMerIn1} and {kMerIn2}.
        @rtype: float
        """

        c = 0.0
        d = 1

        scale1, scale2 = self.__scale(kMerIn1, kMerIn2)

        # Calculate the counter and the denominator of the distance function.
        for i in range(kMerIn1.numberOfKMers) :
            if kMerIn1.kMerCount[i] or kMerIn2.kMerCount[i] :
                c += (abs((scale1 * kMerIn1.kMerCount[i]) -
                          (scale2 * kMerIn2.kMerCount[i])) /
                    ((kMerIn1.kMerCount[i] + 1) * (kMerIn2.kMerCount[i] + 1)))
                d += 1
            #if
        #for

        return c / d
    #__multisetDistance

    def __euclideanDistance(self, kMerIn1, kMerIn2) :
        """
        Calculate the Euclidean distance between two kMer instances.

        @arg kMerIn1: A kMer instance.
        @type kMerIn1: object(kMer)
        @arg kMerIn2: A kMer instance.
        @type kMerIn2: object(kMer)

        @returns: The Euclidean distance between {kMerIn1} and {kMerIn2}.
        @rtype: float
        """

        sumOfSquares = 1

        scale1, scale2 = self.__scale(kMerIn1, kMerIn2)

        # Calculate the counter and the denominator of the distance function.
        for i in range(kMerIn1.numberOfKMers) :
            sumOfSquares += ((scale1 * kMerIn1.kMerCount[i]) -
                (scale2 * kMerIn2.kMerCount[i])) ** 2

        return math.sqrt(sumOfSquares)
    #__euclideanDistance

    def __positiveMultisetDistance(self, kMerIn1, kMerIn2) :
        """
        Calculate the positive multiset distance between two kMer instances.

        @arg kMerIn1: A kMer instance.
        @type kMerIn1: object(kMer)
        @arg kMerIn2: A kMer instance.
        @type kMerIn2: object(kMer)

        @returns: The positive multiset distance between {kMerIn1} and
            {kMerIn2}.
        @rtype: float
        """

        c = 0.0
        d = 1

        scale1, scale2 = self.__positiveScale(kMerIn1, kMerIn2)

        # Calculate the counter and the denominator of the distance function.
        for i in range(kMerIn1.numberOfKMers) :
            if kMerIn1.kMerCount[i] and kMerIn2.kMerCount[i] :
                c += (abs((scale1 * kMerIn1.kMerCount[i]) -
                          (scale2 * kMerIn2.kMerCount[i])) /
                    ((kMerIn1.kMerCount[i] + 1) * (kMerIn2.kMerCount[i] + 1)))
                d += 1
            #if
        #for

        return c / d
    #__positiveMultisetDistance

    def __relativeMultisetDistance(self, kMerIn1, kMerIn2) :
        """
        Calculate the relative multiset distance between two kMer instances.

        @arg kMerIn1: A kMer instance.
        @type kMerIn1: object(kMer)
        @arg kMerIn2: A kMer instance.
        @type kMerIn2: object(kMer)

        @returns: The relative multiset distance between {kMerIn1} and
            {kMerIn2}.
        @rtype: float
        """

        c = 0.0
        d = 1

        scale1, scale2 = self.__scale(kMerIn1, kMerIn2)

        for i in range(kMerIn1.numberOfKMers) :
            for j in range(i) :
                kMerDiff1 = abs(kMerIn1.kMerCount[i] - kMerIn1.kMerCount[j])
                kMerDiff2 = abs(kMerIn2.kMerCount[i] - kMerIn2.kMerCount[j])

                if kMerDiff1 or kMerDiff2 :
                    c += (abs((scale1 * kMerDiff1) - (scale2 * kMerDiff2)) /
                        ((kMerDiff1 + 1) * (kMerDiff2 + 1)))
                    d += 1
                #if
            #for

        return c / d
    #__relativeMultisetDistance
#kMerDiff

def main() :
    """
    Main entry point.
    """

    parser = argparse.ArgumentParser(
        prog = 'kMerDiff',
        formatter_class = argparse.RawDescriptionHelpFormatter,
        description = '',
        epilog = """""")

    parser.add_argument('-i', dest = 'input', type = argparse.FileType('r'),
        required = True, nargs = 2, help = 'The count files.')
    parser.add_argument('-p', dest = 'precision', type = int, default = 3, 
        help = 'Number of decimals.')
    parser.add_argument('-a', dest = 'algorithm', type = int, default = 0,
        help = kMerDiff.algorithmHelp)

    arguments = parser.parse_args()

    kMerDiffInstance = kMerDiff(arguments.algorithm)
    if not kMerDiffInstance.distance :
        print kMerDiff.algorithmError
        parser.print_usage()
        return
    #if

    kMerIn1 = kMer.kMer(0)
    kMerIn2 = kMer.kMer(0)

    kMerIn1.loadKMerCounts(arguments.input[0])
    kMerIn2.loadKMerCounts(arguments.input[1])

    if kMerIn1.kMerLength != kMerIn2.kMerLength :
        print "k-mer lengths of the files differ."
        return
    #if

    print ("%%.%if" % arguments.precision) % \
        kMerDiffInstance.distance(kMerIn1, kMerIn2)
#main

if __name__ == "__main__" :
    main()
