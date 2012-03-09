#!/usr/bin/python

"""
@requires: argparse
@requires: kMer
"""

import argparse
import kMer

def kMerMerge(input1, input2, output) :
    """
    """

    kMerIn = kMer.kMer(0)
    kMerOut = kMer.kMer(0)

    kMerIn.loadKMerCounts(input1)
    kMerOut.loadKMerCounts(input2)

    if kMerIn.kMerLength != kMerOut.kMerLength :
        raise ValueError("k-mer lengths of the files differ.")

    kMerOut.mergeKMerCounts(kMerIn)
    kMerOut.saveKMerCounts(output)
#kMerMerge

def main() :
    """
    Main entry point.
    """

    parser = argparse.ArgumentParser(
        formatter_class = argparse.RawDescriptionHelpFormatter,
        description = '',
        epilog = """""")

    parser.add_argument('-i', dest = 'input', type = argparse.FileType('r'),
        nargs = 2, help = 'The input files.')
    parser.add_argument('-o', dest = 'output', type = argparse.FileType('w'),
        required = True, help = 'The output file name.')

    arguments = parser.parse_args()

    try :
        kMerMerge(arguments.input[0], arguments.input[1], arguments.output)
    except ValueError, error :
        parser.error(error)
#main

if __name__ == "__main__" :
    main()
