#!/usr/bin/python

"""
Toolbox for k-mer profiles.


Use any of the positional arguments with the -h option for more information.
"""

import argparse
import kLib
import kDiffLib

lengthError = "k-mer lengths of the files differ."

def makeProfile(args):
    """
    Make a k-mer profile from a fastq or a fasta file.

    @arg args: Argparse argument list.
    @type args: object
    """
    profile = kLib.kMer()
    profile.analyse(args.input, args.size)
    profile.save(args.output)
#makeProfile

def mergeProfiles(args):
    """
    Merge two k-mer profiles.

    @arg args: Argparse argument list.
    @type args: object
    """
    kMerIn = kLib.kMer()
    kMerOut = kLib.kMer()

    kMerIn.load(args.input1)
    kMerOut.load(args.input2)

    if kMerIn.length != kMerOut.length:
        raise ValueError(lengthError)

    kMerOut.merge(kMerIn)
    kMerOut.save(args.output)
#mergeProfiles

def profileDiff(args):
    """
    Calculate the difference between two k-mer profiles.

    @arg args: Argparse argument list.
    @type args: object
    """
    kDiff = kDiffLib.kMerDiff(args.algorithm, args.function, args.threshold,
        down=args.down)

    profile1 = kLib.kMer(0)
    profile2 = kLib.kMer(0)

    profile1.load(args.input[0])
    profile2.load(args.input[1])

    if profile1.length != profile2.length:
        raise ValueError(lengthError)

    print ("%%.%if" % args.precision) % kDiff.distance(profile1, profile2)
#profileDiff

def diffMatrix(args):
    """
    Make a distance matrix any number of k-mer profiles.

    @arg args: Argparse argument list.
    @type args: object
    """
    if len(args.inputs) < 2:
        raise ValueError("You must give at least two input files.")

    kDiff = kDiffLib.kMerDiff(args.algorithm, args.function, args.threshold,
        down=args.down)

    counts = []
    for i in args.inputs:
        counts.append(kLib.kMer(0))
        counts[-1].load(i)
        if counts[0].length != counts[-1].length:
            raise ValueError(lengthError)
    #for

    kDiffLib.makeDistanceMatrix(counts, args.output, args.precision, kDiff)
#diffMatrix

def smoothProfiles(args):
    """
    Dynamically smooth two k-mer profiles.

    @arg args: Argparse argument list.
    @type args: object
    """
    profile1 = kLib.kMer(0)
    profile2 = kLib.kMer(0)

    profile1.load(args.input[0])
    profile1.load(args.input[1])

    if profile1.length != profile2.length:
        raise ValueError(lengthError)

    dynamicSmooth(profile1, profile2, 0, profile1.number, args.function,
        args.threshold)

    profile1.save(args.output[0])
    profile2.save(args.output[1])
#smoothProfiles

def main():
    """
    Main entry point.
    """
    algorithmHelp = "distance algorithm to use (0=multiset, 1=euclidean, " \
        "2=positive multiset, 3=relative multiset) (default=%(default)s)"
    functionHelp = "smoothing function (default=%(default)s)"
    thresholdHelp = "smoothing threshold (default=%(default)s)"
    inputHelp = "input file name"
    outputHelp = "output file name"
    precisionHelp = "number of decimals (default=%(default)s)"
    downHelp = "scale down"

    parent_parser = argparse.ArgumentParser('parent', add_help=False)

    usage = __doc__.split("\n\n\n")
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=usage[0], epilog=usage[1])
    subparsers = parser.add_subparsers()

    parser_index = subparsers.add_parser("index", parents=[parent_parser],
        description=makeProfile.__doc__.split("\n\n")[0])
    parser_index.add_argument("-i", dest="input", type=argparse.FileType('r'),
        required=True, help=inputHelp)
    parser_index.add_argument("-o", dest="output", type=argparse.FileType('w'),
        required=True, help=outputHelp)
    parser_index.add_argument("-a", dest="inputFormat", default="fastq",
        action="store_const", const="fasta",
        help="use fasta instead of fastq as input format")
    parser_index.add_argument("-r", dest="rc", default=False,
        action="store_true", help="additionally index the reverse complement")
    parser_index.add_argument("-k", dest="size", type=int, required=True,
        help="k-mer size")
    parser_index.set_defaults(func=makeProfile)

    parser_merge = subparsers.add_parser("merge", parents=[parent_parser],
        description=mergeProfiles.__doc__.split("\n\n")[0])
    parser_merge.add_argument("-i", dest="input", type=argparse.FileType('r'),
        nargs=2, help=inputHelp)
    parser_merge.add_argument("-o", dest="output", type=argparse.FileType('w'),
        required=True, help=outputHelp)
    parser_merge.set_defaults(func=mergeProfiles)

    parser_diff = subparsers.add_parser("diff", parents=[parent_parser],
        description=profileDiff.__doc__.split("\n\n")[0])
    parser_diff.add_argument("-i", dest="input", type=argparse.FileType('r'),
        required=True, nargs=2, help="%ss" % inputHelp)
    parser_diff.add_argument("-p", dest="precision", type=int, default=3, 
        help=precisionHelp)
    parser_diff.add_argument("-a", dest="algorithm", type=int, default=0,
        help=algorithmHelp)
    parser_diff.add_argument("-f", dest="function", type=int, default=0,
        help=functionHelp)
    parser_diff.add_argument("-t", dest="threshold", type=int, default=0,
        help=thresholdHelp)
    parser_diff.add_argument("-d", dest="down", default=False,
        action="store_true", help=downHelp)
    parser_diff.set_defaults(func=profileDiff)

    parser_matrix = subparsers.add_parser("matrix", parents=[parent_parser],
        description=diffMatrix.__doc__.split("\n\n")[0])
    parser_matrix.add_argument("-i", dest="inputs", nargs='+',
        type=argparse.FileType('r'), help="%ss" % inputHelp)
    parser_matrix.add_argument("-o", dest="output", required=True,
        type=argparse.FileType('w'), help=outputHelp)
    parser_matrix.add_argument("-p", dest="precision", type=int, default=3,
        help=precisionHelp)
    parser_matrix.add_argument("-a", dest="algorithm", type=int, default=0,
        help=algorithmHelp)
    parser_matrix.add_argument("-f", dest="function", type=int, default=0,
        help=functionHelp)
    parser_matrix.add_argument("-t", dest="threshold", type=int, default=0,
        help=thresholdHelp)
    parser_matrix.add_argument("-s", dest="smooth", default=False,
        action="store_true", help="turn on smoothing")
    parser_matrix.add_argument("-d", dest="down", default=False,
        action="store_true", help=downHelp)
    parser_matrix.set_defaults(func=diffMatrix)

    parser_smooth = subparsers.add_parser("smooth", parents=[parent_parser],
        description=smoothProfiles.__doc__.split("\n\n")[0])
    parser_smooth.add_argument("-i", dest="input", type=argparse.FileType('r'),
        required=True, nargs=2, help="%ss" % inputHelp)
    parser_smooth.add_argument("-o", dest="output", required=True, nargs=2,
        type=argparse.FileType('w'), help="%ss" % outputHelp)
    parser_smooth.add_argument("-f", dest="function", type=int, default=0,
        help=functionHelp)
    parser_smooth.add_argument("-t", dest="threshold", type=int, default=0,
        help=thresholdHelp)
    parser_smooth.set_defaults(func=smoothProfiles)

    arguments = parser.parse_args()

    try:
        arguments.func(arguments)
    except ValueError, error:
        parser.error(error)
#main

if __name__ == "__main__":
    main()
