#!/usr/bin/python

"""
@requires: argparse
@requires: kMer
"""

import argparse
import kMer


def multisetDistance(kMerIn1, kMerIn2) :
    """
    Calculate the multiset distance between two kMer instances.

    @arg kMerIn1: A kMer instance.
    @type kMerIn1: object(kMer)
    @arg kMerIn2: A kMer instance.
    @type kMerIn2: object(kMer)
    """

    scale1 = 1.0
    scale2 = 1.0
    c = 0
    d = 1

    # Calculate the scaling factors in such a way that no element is between 0
    # and 1.
    if kMerIn1.totalKMers < kMerIn2.totalKMers :
        scale1 = float(kMerIn2.totalKMers) / kMerIn1.totalKMers
    else :
        scale2 = float(kMerIn1.totalKMers) / kMerIn2.totalKMers

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
#multisetDistance

def main() :
    """
    Main entry point.
    """

    parser = argparse.ArgumentParser(
        prog = 'kMerDiff',
        formatter_class = argparse.RawDescriptionHelpFormatter,
        description = '',
        epilog = """""")

    parser.add_argument('-i', dest = 'input1', type = argparse.FileType('r'),
        required = True, help = 'The first count file.')
    parser.add_argument('-j', dest = 'input2', type = argparse.FileType('r'),
        required = True, help = 'The first count file.')
    parser.add_argument('-p', dest = 'precision', type = int,
        default = 3, help = 'Number of decimals.')

    arguments = parser.parse_args()

    kMerIn1 = kMer.kMer(0)
    kMerIn2 = kMer.kMer(0)

    kMerIn1.loadKMerCounts(arguments.input1)
    kMerIn2.loadKMerCounts(arguments.input2)

    if kMerIn1.kMerLength != kMerIn2.kMerLength :
        print "k-mer lengths of the files differ."
        return
    #if

    print ("%%.%if" % arguments.precision) % multisetDistance(kMerIn1, kMerIn2)
#main

if __name__ == "__main__" :
    main()
