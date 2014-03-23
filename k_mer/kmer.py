#!/usr/bin/python

"""
Toolbox for k-mer profiles.
"""

import argparse
import sys

import klib
import kdifflib
import metrics

from math import *

from . import docSplit, usage, version

lengthError = "k-mer lengths of the files differ."

def makeProfile(input_handle, output_handle, size):
    """
    Make a k-mer profile from a FASTA file.

    @arg input_handle: Open readable handle to a FASTA file.
    @type input_handle: stream
    @arg output_handle: Open writeable handle to a k-mer profile.
    @type output_handle: stream
    @arg size: Size of k.
    @type size: int
    """
    profile = klib.kMer()
    profile.analyse(input_handle, size)
    profile.save(output_handle)
#makeProfile

def mergeProfiles(input_handles, output_handle):
    """
    Merge two k-mer profiles.

    @arg input_handles: Open readable handles to a pair of k-mer profiles.
    @type input_handles: list(stream)
    @arg output_handle: Open writeable handle to a k-mer profile.
    @type output_handle: stream
    """
    profile1 = klib.kMer()
    profile2 = klib.kMer()

    profile1.load(input_handles[0])
    profile2.load(input_handles[1])

    if profile1.length != profile2.length:
        raise ValueError(lengthError)

    profile2.merge(profile1)
    profile2.save(output_handle)
#mergeProfiles

def balanceProfile(input_handle, output_handle):
    """
    Balance a k-mer profile.

    @arg input_handle: Open readable handle to a k-mer profile.
    @type input_handle: stream
    @arg output_handle: Open writeable handle to a k-mer profile.
    @type output_handle: stream
    """
    profile = klib.kMer()

    profile.load(input_handle)
    profile.balance()
    profile.save(output_handle)
#balanceProfile

def showBalance(input_handle, precision=3):
    """
    Show the balance of a k-mer profile.

    @arg input_handle: Open readable handle to a k-mer profile.
    @type input_handle: stream
    @arg precision: Number of digits in the output.
    @type precision: int
    """
    profile = klib.kMer()

    profile.load(input_handle)
    forward, reverse = profile.split()
    print ("%%.%if" % precision) %  metrics.multisetDistance(forward,
        reverse, metrics.pairwise["diff-prod"])
#showBalance

def showMeanStd(input_handle, precision=3):
    """
    Show the mean and standard deviation of a k-mer profile.

    @arg input_handle: Open readable handle to a k-mer profile.
    @type input_handle: stream
    @arg precision: Number of digits in the output.
    @type precision: int
    """
    profile = klib.kMer()

    profile.load(input_handle)
    print ("%%.%if %%.%if" % (precision, precision)) % \
        metrics.meanStd(profile.count)
#showBalance

def distribution(input_handle, output_handle):
    """
    Calculate the distribution of the values in a k-mer profile.

    @arg input_handle: Open readable handle to a k-mer profile.
    @type input_handle: stream
    @arg output_handle: Open writeable handle to distribution file.
    @type output_handle: stream
    """
    profile = klib.kMer()

    profile.load(input_handle)
    for i in metrics.distribution(profile.count):
        output_handle.write("%i %i\n" % i)
#distribution

def info(input_handle):
    """
    Print some information about the k-mer profile.

    @arg input_handle: Open readable handle to a k-mer profile.
    @type input_handle: stream
    """
    profile = klib.kMer()

    profile.load(input_handle)
    print "k-mer length: %i" % profile.length
    print "Total number of counts: %i" % profile.total
    print "Non-zero counts: %i" % profile.nonZero
#info

def getCount(input_handle, word):
    """
    Retrieve the count for a particular word.

    @arg input_handle: Open readable handle to a k-mer profile.
    @type input_handle: stream
    @arg word: The query word.
    @type word: str
    """
    profile = klib.kMer()

    profile.load(input_handle)
    if profile.length != len(word):
        raise ValueError("The length of the query does not match the profile "
            "length.")
    try:
        offset = profile.DNAToBinary(word)
    except KeyError, error:
        raise ValueError("The input is not a valid DNA sequence.")
    print profile.count[offset]
#getCount

def positiveProfile(input_handles, output_handles):
    """
    Only keep counts that are positive in both profiles.

    @arg input_handles: Open readable handles to a pair of k-mer profiles.
    @type input_handles: list(stream)
    @arg output_handles: Open writeable handles to a pair of k-mer profiles.
    @type output_handles: list(stream)
    """
    profile1 = klib.kMer()
    profile2 = klib.kMer()

    profile1.load(input_handles[0])
    profile2.load(input_handles[1])

    profile1.count = metrics.positive(profile1.count, profile2.count)
    profile2.count = metrics.positive(profile2.count, profile1.count)

    profile1.save(output_handles[0])
    profile2.save(output_handles[1])
#positiveProfile

def scaleProfile(input_handles, output_handles, down=False):
    """
    Scale profiles such that the total number of k-mers is equal.

    @arg input_handles: Open readable handles to a pair of k-mer profiles.
    @type input_handles: list(stream)
    @arg output_handles: Open writeable handles to a pair of k-mer profiles.
    @type output_handles: list(stream)
    @arg down: Scale down.
    @type down: bool
    """
    profile1 = klib.kMer()
    profile2 = klib.kMer()

    profile1.load(input_handles[0])
    profile2.load(input_handles[1])

    scale1, scale2 = metrics.calcScale(profile1.count, profile2.count)
    if down:
        scale1, scale2 = metrics.scaleDown(scale1, scale2)
    metrics.scale(profile1.count, scale1)
    metrics.scale(profile2.count, scale2)

    profile1.save(output_handles[0])
    profile2.save(output_handles[1])
#scaleProfile

def shrinkProfile(input_handle, output_handle, factor):
    """
    Shrink a profile, effectively reducing k.

    @arg input_handle: Open readable handle to a k-mer profile.
    @type input_handle: stream
    @arg output_handle: Open writeable handle to a k-mer profile.
    @type output_handle: stream
    @arg factor: Scaling factor.
    @type factor: int
    """
    profile = klib.kMer()

    profile.load(input_handle)
    profile.shrink(factor)
    profile.save(output_handle)
#shrinkProfile

def shuffleProfile(input_handle, output_handle):
    """
    Randomise a profile.

    @arg input_handle: Open readable handle to a k-mer profile.
    @type input_handle: stream
    @arg output_handle: Open writeable handle to a k-mer profile.
    @type output_handle: stream
    """
    profile = klib.kMer()

    profile.load(input_handle)
    profile.shuffle()
    profile.save(output_handle)
#shuffleProfile

def smoothProfile(input_handles, output_handles, summary, summary_func="",
        threshold=0):
    """
    Smooth two profiles by collapsing sub-profiles.

    @arg input_handles: Open readable handles to a pair of k-mer profiles.
    @type input_handles: list(stream)
    @arg output_handles: Open writeable handles to a pair of k-mer profiles.
    @type output_handles: list(stream)
    @arg summary: Name of the summary function.
    @type summary: str
    @arg summary_func: Custom summary function.
    @type summary_func: str
    @arg threshold: Threshold for the summary function.
    @type threshold: int
    """
    smooth_function = metrics.summary[summary]
    if summary_func:
        smooth_function = lambda x: eval(summary_func)

    diff = kdifflib.kMerDiff(summary=smooth_function, threshold=threshold)
    profile1 = klib.kMer()
    profile2 = klib.kMer()

    profile1.load(input_handles[0])
    profile2.load(input_handles[1])

    diff.dynamicSmooth(profile1, profile2)

    profile1.save(output_handles[0])
    profile2.save(output_handles[1])
#smoothProfile

def diffProfile(input_handles, euclidean=False, pairwise="diff-prod",
        pairwise_func="", smooth=False, summary="min", summary_func="",
        threshold=0, scale=False, down=False, positive=False, balance=False,
        precision=3):
    """
    Calculate the difference between two k-mer profiles.

    @arg input_handles: Open readable handles to a pair of k-mer profiles.
    @type input_handles: list(stream)
    @arg euclidean: Use the Euclidean distance.
    @type euclidean: bool
    @arg pairwise: Name of the pairwise distance function.
    @type pairwise: str
    @arg pairwise_func: Custom pairwise distance function.
    @type pairwise_func: str
    @arg smooth: Enable smoothing.
    @type smooth: bool
    @arg summary: Name of the summary function.
    @type summary: str
    @arg summary_func: Custom summary function.
    @type summary_func: str
    @arg threshold: Threshold for the summary function.
    @type threshold: int
    @arg scale: Scale the profiles.
    @type scale: bool
    @arg down: Scale down.
    @type down: bool
    @arg positive: Only use positive values.
    @type positive: bool
    @arg balance: Balance the profiles.
    @type balance: bool
    @arg precision: Number of digits in the output.
    @type precision: int
    """
    summary_function = metrics.summary[summary]
    if summary_func:
        summary_function = lambda x: eval(summary_func)

    pairwise_function = metrics.pairwise[pairwise]
    if pairwise_func:
        pairwise_function = lambda x, y: eval(pairwise_func)

    diff = kdifflib.kMerDiff(balance=balance, positive=positive, smooth=smooth,
        summary=summary_function, threshold=threshold, scale=scale,
        scaleDown=down, multiset=not euclidean, pairwise=pairwise_function)

    profile1 = klib.kMer()
    profile2 = klib.kMer()

    profile1.load(input_handles[0])
    profile2.load(input_handles[1])

    if profile1.length != profile2.length:
        raise ValueError(lengthError)

    print ("%%.%if" % precision) % diff.calcDistance(profile1, profile2)
#diffProfile

def diffMatrix(input_handles, euclidean=False, pairwise="diff-prod",
        pairwise_func="", smooth=False, summary="min", summary_func="",
        threshold=0, scale=False, down=False, positive=False, balance=False,
        precision=3):
    """
    Make a distance matrix any number of k-mer profiles.

    @arg input_handles: Open readable handles to a list of k-mer profiles.
    @type input_handles: list(stream)
    @arg euclidean: Use the Euclidean distance.
    @type euclidean: bool
    @arg pairwise: Name of the pairwise distance function.
    @type pairwise: str
    @arg pairwise_func: Custom pairwise distance function.
    @type pairwise_func: str
    @arg smooth: Enable smoothing.
    @type smooth: bool
    @arg summary: Name of the summary function.
    @type summary: str
    @arg summary_func: Custom summary function.
    @type summary_func: str
    @arg threshold: Threshold for the summary function.
    @type threshold: int
    @arg scale: Scale the profiles.
    @type scale: bool
    @arg down: Scale down.
    @type down: bool
    @arg positive: Only use positive values.
    @type positive: bool
    @arg balance: Balance the profiles.
    @type balance: bool
    @arg precision: Number of digits in the output.
    @type precision: int
    """
    if len(input_handles) < 2:
        raise ValueError("You must give at least two input files.")

    summary_function = metrics.summary[summary]
    if summary_func:
        summary_function = lambda x: eval(summary_func)

    pairwise_function = metrics.pairwise[pairwise]
    if pairwise_func:
        pairwise_function = lambda x, y: eval(pairwise_func)

    diff = kdifflib.kMerDiff(balance=balance, positive=positive,
        smooth=smooth, summary=summary_function, threshold=threshold,
        scale=scale, scaleDown=down, multiset=not euclidean,
        pairwise=pairwise_function)

    counts = []
    for i in input_handles:
        counts.append(klib.kMer())
        counts[-1].load(i)
        if counts[0].length != counts[-1].length:
            raise ValueError(lengthError)
    #for

    kdifflib.makeDistanceMatrix(counts, output_handle, precision, diff)
#diffMatrix

def main():
    """
    Main entry point.
    """
    output_parser = argparse.ArgumentParser(add_help=False)
    output_parser.add_argument("output_handle", metavar="OUTPUT",
        type=argparse.FileType('w'), help="output file")

    input_parser = argparse.ArgumentParser(add_help=False)
    input_parser.add_argument("input_handle", metavar="INPUT",
        type=argparse.FileType('r'), help="input file")

    pairOut_parser = argparse.ArgumentParser(add_help=False)
    pairOut_parser.add_argument("output_handles", metavar="OUTPUT", nargs=2,
        type=argparse.FileType('w'), help="pair of output files")

    pairIn_parser = argparse.ArgumentParser(add_help=False)
    pairIn_parser.add_argument("input_handles", metavar="INPUT", nargs=2,
        type=argparse.FileType('r'), help="pair of input files")

    listIn_parser = argparse.ArgumentParser(add_help=False)
    listIn_parser.add_argument("input_handles", metavar="INPUT", nargs='+',
        type=argparse.FileType('r'), help="list of input files")

    scale_parser = argparse.ArgumentParser(add_help=False)
    scale_parser.add_argument("-d", dest="down", default=False,
        action="store_true", help="scale down (default=%(default)s)")

    smooth_parser = argparse.ArgumentParser(add_help=False)
    smooth_parser.add_argument("-s", dest="summary", type=str, default="min",
        choices=metrics.summary, help="summary function for dynamic smoothing "
        '(%(type)s default="%(default)s")')
    smooth_parser.add_argument("--summary-function", dest="summary_func",
        type=str, default="", help="custom summary function "
        '(%(type)s default="%(default)s")')
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
    diff_parser.add_argument("-P", dest="pairwise", type=str,
        default="diff-prod", choices=metrics.pairwise,
        help="paiwise distance function for the multiset distance "
        '(%(type)s default="%(default)s")')
    diff_parser.add_argument("--pairwise-function", dest="pairwise_func",
        type=str, default="", help="custom pairwise function "
        '(%(type)s default="%(default)s")')

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=usage[0], epilog=usage[1])
    parser.add_argument('-v', action="version", version=version(parser.prog))
    subparsers = parser.add_subparsers()

    parser_index = subparsers.add_parser("index", parents=[input_parser,
        output_parser], description=docSplit(makeProfile))
    parser_index.add_argument("size", metavar="SIZE", type=int,
        help="k-mer size (%(type)s)")
    parser_index.set_defaults(func=makeProfile)

    parser_merge = subparsers.add_parser("merge", parents=[pairIn_parser,
        output_parser], description=docSplit(mergeProfiles))
    parser_merge.set_defaults(func=mergeProfiles)

    parser_balance = subparsers.add_parser("balance", parents=[input_parser,
        output_parser], description=docSplit(balanceProfile))
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
        parents=[input_parser, output_parser],
        description=docSplit(distribution))
    parser_distr.set_defaults(func=distribution)

    parser_info = subparsers.add_parser("info", 
        parents=[input_parser],
        description=docSplit(info))
    parser_info.set_defaults(func=info)

    parser_getcount = subparsers.add_parser("getcount", 
        parents=[input_parser],
        description=docSplit(getCount))
    parser_getcount.add_argument("word", metavar="WORD", type=str,
        help="the word in question (%(type)s)")
    parser_getcount.set_defaults(func=getCount)

    parser_positive = subparsers.add_parser("positive",
        parents=[pairIn_parser, pairOut_parser],
        description=docSplit(positiveProfile))
    parser_positive.set_defaults(func=positiveProfile)

    parser_scale = subparsers.add_parser("scale", parents=[pairIn_parser,
        pairOut_parser, scale_parser], description=docSplit(scaleProfile))
    parser_scale.set_defaults(func=scaleProfile)

    parser_shrink = subparsers.add_parser("shrink", parents=[input_parser,
        output_parser], description=docSplit(shrinkProfile))
    parser_shrink.add_argument("factor", metavar="FACTOR", type=int,
        help="shrinking factor")
    parser_shrink.set_defaults(func=shrinkProfile)

    parser_shrink = subparsers.add_parser("shuffle", parents=[input_parser,
        output_parser], description=docSplit(shuffleProfile))
    parser_shrink.set_defaults(func=shuffleProfile)

    parser_smooth = subparsers.add_parser("smooth", parents=[pairIn_parser,
        pairOut_parser, smooth_parser], description=docSplit(smoothProfile))
    parser_smooth.set_defaults(func=smoothProfile)

    parser_diff = subparsers.add_parser("diff", parents=[diff_parser,
        pairIn_parser], description=docSplit(diffProfile))
    parser_diff.set_defaults(func=diffProfile)

    parser_matrix = subparsers.add_parser("matrix", parents=[diff_parser,
        listIn_parser, output_parser], description=docSplit(diffMatrix))
    parser_matrix.set_defaults(func=diffMatrix)

    arguments = parser.parse_args()

    try:
        arguments.func(**{k: v for k, v in vars(arguments).items()
            if k not in ("func", "subcommand")})
    except ValueError, error:
        parser.error(error)
#main

if __name__ == "__main__":
    main()
