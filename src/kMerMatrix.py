#!/usr/bin/python

"""
@requires: argparse
@requires: kMer
@requires: kMerDiff
"""

import argparse
import kMer
import kMerDiff


def makeDistanceMatrix(kMerCounts, output, precision, kMerDiffInstance) :
    """
    Make a distance matrix for all the kMer instances.

    @arg kMerCounts: List of kMer instances.
    @type kMerCounts: list(kMer)
    @arg output: Open handle to a writable file.
    @type output: stream
    @arg precision: Number of decimals.
    @type precision: integer
    @arg kMerDiffInstance: A kMerDiff object.
    @type kMerDiffInstance: list(kMerDiff)
    """

    numberOfInputs = len(kMerCounts)
    output.write("%i\n" % numberOfInputs)
    for i in range(1, numberOfInputs + 1) :
        output.write("%i\n" % i)
    for i in range(1, numberOfInputs) :
        for j in range(i) :
            if (j) :
                output.write(' ')
            output.write(("%%.%if" % precision) % 
                kMerDiffInstance.distance(kMerCounts[i], kMerCounts[j]))
        #for
        output.write('\n')
    #for
#makeDistanceMatrix

def main() :
    """
    Main entry point.
    """

    parser = argparse.ArgumentParser(
        prog = 'kMerMatrix',
        formatter_class = argparse.RawDescriptionHelpFormatter,
        description = '',
        epilog = """""")

    parser.add_argument('-i', dest = 'input', type = argparse.FileType('r'),
        nargs = '+', help = 'The first count file.')
    parser.add_argument('-o', dest = 'output', type = argparse.FileType('w'),
        required = True, help = 'The output file name.')
    parser.add_argument('-p', dest = 'precision', type = int,
        default = 3, help = 'Number of decimals.')
    parser.add_argument('-a', dest = 'algorithm', type = int, default = 0,
        help = kMerDiff.kMerDiff.algorithmHelp)

    arguments = parser.parse_args()

    if len(arguments.input) < 2 :
        print "You must give at least two input files."
        return
    #if

    kMerDiffInstance = kMerDiff.kMerDiff(arguments.algorithm)
    if not kMerDiffInstance.distance :
        print kMerDiff.kMerDiff.algorithmError
        parser.print_usage()
        return
    #if

    kMerCounts = []
    for i in arguments.input :
        kMerCounts.append(kMer.kMer(0))
        kMerCounts[-1].loadKMerCounts(i)
        if kMerCounts[0].kMerLength != kMerCounts[-1].kMerLength :
            print "k-mer lengths of the files differ."
            return
        #if
    #for

    makeDistanceMatrix(kMerCounts, arguments.output, arguments.precision,
        kMerDiffInstance)
#main

if __name__ == "__main__" :
    main()
