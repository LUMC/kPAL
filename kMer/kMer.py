#!/usr/bin/python

"""
Toolbox for k-mer profiles.


Use any of the positional arguments with the -h option for more information.
"""

import argparse
import sys

import kLib
import kDiffLib
import metrics

lengthError = "k-mer lengths of the files differ."

def docSplit(func):
    """
    Return the header of the docstring of a function.

    @arg func: A function.
    @type func: function
    """
    return func.__doc__.split("\n\n")[0]
#docSplit

def makeProfile(args):
    """
    Make a k-mer profile from a fasta file.

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
    profile1 = kLib.kMer()
    profile2 = kLib.kMer()

    profile1.load(args.input[0])
    profile2.load(args.input[1])

    if profile1.length != profile2.length:
        raise ValueError(lengthError)

    profile2.merge(profile1)
    profile2.save(args.output)
#mergeProfiles

def balanceProfile(args):
    """
    Balance a k-mer profile.

    @arg args: Argparse argument list.
    @type args: object
    """
    profile = kLib.kMer()

    profile.load(args.input)
    profile.balance()
    profile.save(args.output)
#balanceProfile

def showBalance(args):
    """
    Show the balance of a k-mer profile.

    @arg args: Argparse argument list.
    @type args: object
    """
    profile = kLib.kMer()

    profile.load(args.input)
    forward, reverse = profile.split()
    print ("%%.%if" % args.precision) %  metrics.multisetDistance(forward,
        reverse, metrics.pairwise[0])
#showBalance

def showMeanStd(args):
    """
    Show the mean and standard deviation of a k-mer profile.

    @arg args: Argparse argument list.
    @type args: object
    """
    profile = kLib.kMer()

    profile.load(args.input)
    print ("%%.%if %%.%if" % (args.precision, args.precision)) % \
        metrics.meanStd(profile.count)
#showBalance

def distribution(args):
    """
    Calculate the distribution of the values in a k-mer profile.

    @arg args: Argparse argument list.
    @type args: object
    """
    profile = kLib.kMer()

    profile.load(args.input)
    for i in metrics.distribution(profile.count):
        args.output.write("%i %i\n" % i)
#distribution

def info(args):
    """
    Print some information about the k-mer profile.

    @arg args: Argparse argument list.
    @type args: object
    """
    profile = kLib.kMer()

    profile.load(args.input)
    print "k-mer length: %i" % profile.length
    print "Total number of counts: %i" % profile.total
    print "Non-zero counts: %i" % profile.nonZero
#info

def getCount(args):
    """
    Retrieve the count for a particular word.

    @arg args: Argparse argument list.
    @type args: object
    """
    profile = kLib.kMer()

    profile.load(args.input)
    if profile.length != len(args.word):
        raise ValueError("The length of the query does not match the profile "
            "length.")
    try:
        offset = profile.DNAToBinary(args.word)
    except KeyError, error:
        raise ValueError("The input is not a valid DNA sequence.")
    print profile.count[offset]
#getCount

def positiveProfile(args):
    """
    Only keep counts that are positive in both profiles.

    @arg args: Argparse argument list.
    @type args: object
    """
    profile1 = kLib.kMer()
    profile2 = kLib.kMer()

    profile1.load(args.input[0])
    profile2.load(args.input[1])

    profile1.count = metrics.positive(profile1.count, profile2.count)
    profile2.count = metrics.positive(profile2.count, profile1.count)

    profile1.save(args.output[0])
    profile2.save(args.output[1])
#positiveProfile

def scaleProfile(args):
    """
    Scale profiles such that the total number of k-mers is equal.

    @arg args: Argparse argument list.
    @type args: object
    """
    profile1 = kLib.kMer()
    profile2 = kLib.kMer()

    profile1.load(args.input[0])
    profile2.load(args.input[1])

    scale1, scale2 = metrics.calcScale(profile1.count, profile2.count)
    if args.down:
        scale1, scale2 = metrics.scaleDown(scale1, scale2)
    metrics.scale(profile1.count, scale1)
    metrics.scale(profile2.count, scale2)

    profile1.save(args.output[0])
    profile2.save(args.output[1])
#scaleProfile

def shrinkProfile(args):
    """
    Shrink a profile, effectively reducing k.

    @arg args: Argparse argument list.
    @type args: object
    """
    profile = kLib.kMer()

    profile.load(args.input)
    profile.shrink(args.factor)
    profile.save(args.output)
#shrinkProfile

def shuffleProfile(args):
    """
    Randomise a profile.

    @arg args: Argparse argument list.
    @type args: object
    """
    profile = kLib.kMer()

    profile.load(args.input)
    profile.shuffle()
    profile.save(args.output)
#shuffleProfile

def smoothProfile(args):
    """
    Smooth two profiles by collapsing sub-profiles.

    @arg args: Argparse argument list.
    @type args: object
    """
    diff = kDiffLib.kMerDiff(summary=args.summary, threshold=args.threshold)
    profile1 = kLib.kMer()
    profile2 = kLib.kMer()

    profile1.load(args.input[0])
    profile2.load(args.input[1])

    diff.dynamicSmooth(profile1, profile2)

    profile1.save(args.output[0])
    profile2.save(args.output[1])
#smoothProfile

def diffProfile(args):
    """
    Calculate the difference between two k-mer profiles.

    @arg args: Argparse argument list.
    @type args: object
    """
    diff = kDiffLib.kMerDiff(balance=args.balance, positive=args.positive,
        smooth=args.smooth, summary=args.summary, threshold=args.threshold,
        scale=args.scale, scaleDown=args.down, multiset=not args.euclidean,
        pairwise=args.pairwise)

    profile1 = kLib.kMer()
    profile2 = kLib.kMer()

    profile1.load(args.input[0])
    profile2.load(args.input[1])

    if profile1.length != profile2.length:
        raise ValueError(lengthError)

    print ("%%.%if" % args.precision) % diff.calcDistance(profile1, profile2)
#diffProfile

def diffMatrix(args):
    """
    Make a distance matrix any number of k-mer profiles.

    @arg args: Argparse argument list.
    @type args: object
    """
    if len(args.inputs) < 2:
        raise ValueError("You must give at least two input files.")

    diff = kDiffLib.kMerDiff(balance=args.balance, positive=args.positive,
        smooth=args.smooth, summary=args.summary, threshold=args.threshold,
        scale=args.scale, scaleDown=args.down, multiset=not args.euclidean,
        pairwise=args.pairwise)

    counts = []
    for i in args.inputs:
        counts.append(kLib.kMer())
        counts[-1].load(i)
        if counts[0].length != counts[-1].length:
            raise ValueError(lengthError)
    #for

    kDiffLib.makeDistanceMatrix(counts, args.output, args.precision, diff)
#diffMatrix

def main():
    """
    Main entry point.
    """
    output_parser = argparse.ArgumentParser(add_help=False)
    output_parser.add_argument("-o", dest="output", default=sys.stdout,
        type=argparse.FileType('w'), help="output file (default=<stdout>)")

    input_parser = argparse.ArgumentParser(add_help=False)
    input_parser.add_argument("-i", dest="input", default=sys.stdin,
        type=argparse.FileType('r'), help="input file (default=<stdin>)")

    pairOut_parser = argparse.ArgumentParser(add_help=False)
    pairOut_parser.add_argument("-o", dest="output", nargs=2, required=True,
        type=argparse.FileType('w'), help="pair of output files")

    pairIn_parser = argparse.ArgumentParser(add_help=False)
    pairIn_parser.add_argument("-i", dest="input", nargs=2, required=True, 
        type=argparse.FileType('r'), help="pair of input files")

    scale_parser = argparse.ArgumentParser(add_help=False)
    scale_parser.add_argument("-d", dest="down", default=False,
        action="store_true", help="scale down (default=%(default)s)")

    smooth_parser = argparse.ArgumentParser(add_help=False)
    smooth_parser.add_argument("-s", dest="summary", type=int, default=0,
        help="summary function for dynamic smoothing "
        "(%(type)s default=%(default)s)")
    smooth_parser.add_argument("-t", dest="threshold", type=int, default=0,
        help="threshold for the summary function "
        "(%(type)s default=%(default)s)")

    precision_parser = argparse.ArgumentParser(add_help=False)
    precision_parser.add_argument("-n", dest="precision", type=int, default=3,
        help="number of decimals (%(type)s default=%(default)s)")

    diff_parser = argparse.ArgumentParser(add_help=False,
        parents=[scale_parser, smooth_parser, precision_parser])
    diff_parser.add_argument("-b", dest="balance", default=False,
        action="store_true", help="balance the profiles (default=%(default)s)")
    diff_parser.add_argument("-p", dest="positive", default=False,
        action="store_true", help="use only positive values "
        "(default=%(default)s)")
    diff_parser.add_argument("-S", dest="scale", default=False,
        action="store_true", help="scale the profiles (default=%(default)s)")
    diff_parser.add_argument("-m", dest="smooth", default=False,
        action="store_true", help="smooth the profiles (default=%(default)s)")
    diff_parser.add_argument("-e", dest="euclidean", default=False,
        action="store_true", help="use the euclidean distance metric "
        "(default=%(default)s)")
    diff_parser.add_argument("-P", dest="pairwise", type=int, default=0,
        help="paiwise distance function for the multiset distance "
        "(%(type)s default=%(default)s)")

    usage = __doc__.split("\n\n\n")
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=usage[0], epilog=usage[1])
    subparsers = parser.add_subparsers()

    parser_index = subparsers.add_parser("index", parents=[output_parser],
        description=docSplit(makeProfile))
    parser_index.add_argument("-i", dest="input", type=argparse.FileType('r'),
        default=sys.stdin, help="input file in fasta format (default=<stdin>)")
    parser_index.add_argument("-k", dest="size", type=int, required=True,
        help="k-mer size (%(type)s)")
    parser_index.set_defaults(func=makeProfile)

    parser_merge = subparsers.add_parser("merge", parents=[output_parser,
        pairIn_parser], description=docSplit(mergeProfiles))
    parser_merge.set_defaults(func=mergeProfiles)

    parser_balance = subparsers.add_parser("balance", parents=[output_parser,
        input_parser], description=docSplit(balanceProfile))
    parser_balance.set_defaults(func=balanceProfile)

    parser_balance = subparsers.add_parser("showbalance", 
        parents=[input_parser, precision_parser],
        description=docSplit(showBalance))
    parser_balance.set_defaults(func=showBalance)

    parser_meanstd = subparsers.add_parser("meanstd", 
        parents=[input_parser, precision_parser],
        description=docSplit(showMeanStd))
    parser_meanstd.set_defaults(func=showMeanStd)

    parser_distr = subparsers.add_parser("distr", 
        parents=[output_parser, input_parser],
        description=docSplit(distribution))
    parser_distr.set_defaults(func=distribution)

    parser_info = subparsers.add_parser("info", 
        parents=[input_parser],
        description=docSplit(info))
    parser_info.set_defaults(func=info)

    parser_getcount = subparsers.add_parser("getcount", 
        parents=[input_parser],
        description=docSplit(getCount))
    parser_getcount.add_argument("-w", dest="word", type=str, required=True,
        help="the word in question (%(type)s)")
    parser_getcount.set_defaults(func=getCount)

    parser_positive = subparsers.add_parser("positive",
        parents=[pairOut_parser, pairIn_parser],
        description=docSplit(positiveProfile))
    parser_positive.set_defaults(func=positiveProfile)

    parser_scale = subparsers.add_parser("scale", parents=[pairOut_parser,
        pairIn_parser, scale_parser], description=docSplit(scaleProfile))
    parser_scale.set_defaults(func=scaleProfile)

    parser_shrink = subparsers.add_parser("shrink", parents=[output_parser,
        input_parser], description=docSplit(shrinkProfile))
    parser_shrink.add_argument("-f", dest="factor", type=int, required=True,
        help="shrinking factor")
    parser_shrink.set_defaults(func=shrinkProfile)

    parser_shrink = subparsers.add_parser("shuffle", parents=[output_parser,
        input_parser], description=docSplit(shuffleProfile))
    parser_shrink.set_defaults(func=shuffleProfile)

    parser_smooth = subparsers.add_parser("smooth", parents=[pairOut_parser,
        pairIn_parser, smooth_parser], description=docSplit(smoothProfile))
    parser_smooth.set_defaults(func=smoothProfile)

    parser_diff = subparsers.add_parser("diff", parents=[pairIn_parser,
        diff_parser], description=docSplit(diffProfile))
    parser_diff.set_defaults(func=diffProfile)

    parser_matrix = subparsers.add_parser("matrix", parents=[diff_parser,
        output_parser], description=docSplit(diffMatrix))
    parser_matrix.add_argument("-i", dest="inputs", nargs='+', required=True,
        type=argparse.FileType('r'), help="list of input files")
    parser_matrix.set_defaults(func=diffMatrix)

    arguments = parser.parse_args()

    try:
        arguments.func(arguments)
    except ValueError, error:
        parser.error(error)
#main

if __name__ == "__main__":
    main()
