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

def kMerMatrix(inputs, output, algorithm, precision, down) :
    """
    """

    if len(inputs) < 2 :
        raise ValueError("You must give at least two input files.")

    kMerDiffInstance = kMerDiff.kMerDiff(algorithm, down=down)
    if not kMerDiffInstance.distance :
        raise ValueError(kMerDiff.kMerDiff.algorithmError)

    kMerCounts = []
    for i in inputs :
        kMerCounts.append(kMer.kMer(0))
        kMerCounts[-1].loadKMerCounts(i)
        if kMerCounts[0].kMerLength != kMerCounts[-1].kMerLength :
            raise ValueError("k-mer lengths of the files differ.")
    #for

    makeDistanceMatrix(kMerCounts, output, precision, kMerDiffInstance)
#kMerMatrix

def main() :
    """
    Main entry point.
    """

    parser = argparse.ArgumentParser(
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
    parser.add_argument('-d', dest = 'down', default = False,
        action = 'store_true', help = kMerDiff.kMerDiff.downHelp)

    arguments = parser.parse_args()

    try :
        kMerMatrix(arguments.input, arguments.output, arguments.algorithm,
            arguments.precision, arguments.down)
    except ValueError, error :
        parser.error(error)
#main

if __name__ == "__main__" :
    main()
