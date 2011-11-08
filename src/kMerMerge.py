#!/usr/bin/python

"""
@requires: argparse
@requires: kMer
"""

import argparse
import kMer

def main() :
    """
    Main entry point.
    """

    parser = argparse.ArgumentParser(
        prog = 'kMerMerge',
        formatter_class = argparse.RawDescriptionHelpFormatter,
        description = '',
        epilog = """""")

    parser.add_argument('-i', dest = 'input1', type = argparse.FileType('r'),
        required = True, help = 'The first count file.')
    parser.add_argument('-j', dest = 'input2', type = argparse.FileType('r'),
        required = True, help = 'The first count file.')
    parser.add_argument('-o', dest = 'output', type = argparse.FileType('w'),
        required = True, help = 'The output file name.')

    arguments = parser.parse_args()

    kMerIn = kMer.kMer(0)
    kMerOut = kMer.kMer(0)

    kMerIn.loadKMerCounts(arguments.input1)
    kMerOut.loadKMerCounts(arguments.input2)

    if kMerIn.kMerLength == kMerOut.kMerLength :
        kMerOut.mergeKMerCounts(kMerIn)
        kMerOut.saveKMerCounts(arguments.output)
    #if
    else :
        print "k-mer lengths of the files differ."
#main

if __name__ == "__main__" :
    main()
