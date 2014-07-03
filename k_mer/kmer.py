#!/usr/bin/python

"""
Toolbox for k-mer profiles.
"""

from __future__ import division

import argparse
import sys

from math import *

from . import (ProtectedFileType, doc_split, usage, version, klib, kdifflib,
    metrics)

length_error = "k-mer lengths of the files differ."

def index(input_handle, output_handle, size):
    """
    Make a k-mer profile from a FASTA file.

    :arg input_handle: Open readable handle to a FASTA file.
    :type input_handle: stream
    :arg output_handle: Open writeable handle to a k-mer profile.
    :type output_handle: stream
    :arg size: Size of k.
    :type size: int
    """
    profile = klib.kMer()
    profile.analyse(input_handle, size)
    profile.save(output_handle)
#index

def merge(input_handles, output_handle, merger="sum", merge_func=""):
    """
    Merge two k-mer profiles.

    :arg input_handles: Open readable handles to a pair of k-mer profiles.
    :type input_handles: list(stream)
    :arg output_handle: Open writeable handle to a k-mer profile.
    :type output_handle: stream
    :arg merger: Merge function.
    :type merger: function
    :arg merge_func: Custom merge function.
    :type merge_func: str
    """
    merge_function = metrics.mergers[merger]
    if merge_func:
        merge_function = eval("lambda " + merge_func)

    profile1 = klib.kMer()
    profile2 = klib.kMer()

    profile1.load(input_handles[0])
    profile2.load(input_handles[1])

    if profile1.length != profile2.length:
        raise ValueError(length_error)

    profile2.merge(profile1, merge_function)
    profile2.save(output_handle)
#merge

def balance(input_handle, output_handle):
    """
    Balance a k-mer profile.

    :arg input_handle: Open readable handle to a k-mer profile.
    :type input_handle: stream
    :arg output_handle: Open writeable handle to a k-mer profile.
    :type output_handle: stream
    """
    profile = klib.kMer()

    profile.load(input_handle)
    profile.balance()
    profile.save(output_handle)
#balance

def get_balance(input_handle, precision=3):
    """
    Show the balance of a k-mer profile.

    :arg input_handle: Open readable handle to a k-mer profile.
    :type input_handle: stream
    :arg precision: Number of digits in the output.
    :type precision: int
    """
    profile = klib.kMer()

    profile.load(input_handle)
    forward, reverse = profile.split()
    print ("%%.%if" % precision) %  metrics.multiset(forward, reverse,
        metrics.pairwise["diff-prod"])
#get_balance

def get_stats(input_handle, precision=3):
    """
    Show the mean and standard deviation of a k-mer profile.

    :arg input_handle: Open readable handle to a k-mer profile.
    :type input_handle: stream
    :arg precision: Number of digits in the output.
    :type precision: int
    """
    profile = klib.kMer()

    profile.load(input_handle)
    print ("%%.%if %%.%if" % (precision, precision)) % \
        metrics.stats(profile.count)
#get_stats

def distribution(input_handle, output_handle):
    """
    Calculate the distribution of the values in a k-mer profile.

    :arg input_handle: Open readable handle to a k-mer profile.
    :type input_handle: stream
    :arg output_handle: Open writeable handle to distribution file.
    :type output_handle: stream
    """
    profile = klib.kMer()

    profile.load(input_handle)
    for i in metrics.distribution(profile.count):
        output_handle.write("%i %i\n" % i)
#distribution

def info(input_handle):
    """
    Print some information about the k-mer profile.

    :arg input_handle: Open readable handle to a k-mer profile.
    :type input_handle: stream
    """
    profile = klib.kMer()

    profile.load(input_handle)
    print "k-mer length: %i" % profile.length
    print "Total number of counts: %i" % profile.total
    print "Non-zero counts: %i" % profile.non_zero
#info

def get_count(input_handle, word):
    """
    Retrieve the count for a particular word.

    :arg input_handle: Open readable handle to a k-mer profile.
    :type input_handle: stream
    :arg word: The query word.
    :type word: str
    """
    profile = klib.kMer()

    profile.load(input_handle)
    if profile.length != len(word):
        raise ValueError("The length of the query does not match the profile "
            "length.")
    try:
        offset = profile.dna_to_binary(word)
    except KeyError, error:
        raise ValueError("The input is not a valid DNA sequence.")
    print profile.count[offset]
#get_count

def positive(input_handles, output_handles):
    """
    Only keep counts that are positive in both profiles.

    :arg input_handles: Open readable handles to a pair of k-mer profiles.
    :type input_handles: list(stream)
    :arg output_handles: Open writeable handles to a pair of k-mer profiles.
    :type output_handles: list(stream)
    """
    profile1 = klib.kMer()
    profile2 = klib.kMer()

    profile1.load(input_handles[0])
    profile2.load(input_handles[1])

    profile1.count = metrics.positive(profile1.count, profile2.count)
    profile2.count = metrics.positive(profile2.count, profile1.count)

    profile1.save(output_handles[0])
    profile2.save(output_handles[1])
#positive

def scale(input_handles, output_handles, down=False):
    """
    Scale profiles such that the total number of k-mers is equal.

    :arg input_handles: Open readable handles to a pair of k-mer profiles.
    :type input_handles: list(stream)
    :arg output_handles: Open writeable handles to a pair of k-mer profiles.
    :type output_handles: list(stream)
    :arg down: Scale down.
    :type down: bool
    """
    profile1 = klib.kMer()
    profile2 = klib.kMer()

    profile1.load(input_handles[0])
    profile2.load(input_handles[1])

    scale1, scale2 = metrics.get_scale(profile1.count, profile2.count)
    if down:
        scale1, scale2 = metrics.scale_down(scale1, scale2)
    profile1.count = metrics.scale(profile1.count, scale1)
    profile2.count = metrics.scale(profile2.count, scale2)

    profile1.save(output_handles[0])
    profile2.save(output_handles[1])
#scale

def shrink(input_handle, output_handle, factor):
    """
    Shrink a profile, effectively reducing k.

    :arg input_handle: Open readable handle to a k-mer profile.
    :type input_handle: stream
    :arg output_handle: Open writeable handle to a k-mer profile.
    :type output_handle: stream
    :arg factor: Scaling factor.
    :type factor: int
    """
    profile = klib.kMer()

    profile.load(input_handle)
    profile.shrink(factor)
    profile.save(output_handle)
#shrink

def shuffle(input_handle, output_handle):
    """
    Randomise a profile.

    :arg input_handle: Open readable handle to a k-mer profile.
    :type input_handle: stream
    :arg output_handle: Open writeable handle to a k-mer profile.
    :type output_handle: stream
    """
    profile = klib.kMer()

    profile.load(input_handle)
    profile.shuffle()
    profile.save(output_handle)
#shuffle

def smooth(input_handles, output_handles, summary, summary_func="",
        threshold=0):
    """
    Smooth two profiles by collapsing sub-profiles.

    :arg input_handles: Open readable handles to a pair of k-mer profiles.
    :type input_handles: list(stream)
    :arg output_handles: Open writeable handles to a pair of k-mer profiles.
    :type output_handles: list(stream)
    :arg summary: Name of the summary function.
    :type summary: str
    :arg summary_func: Custom summary function.
    :type summary_func: str
    :arg threshold: Threshold for the summary function.
    :type threshold: int
    """
    smooth_function = metrics.summary[summary]
    if summary_func:
        smooth_function = eval("lambda " + summary_func)

    diff = kdifflib.kMerDiff(summary=smooth_function, threshold=threshold)
    profile1 = klib.kMer()
    profile2 = klib.kMer()

    profile1.load(input_handles[0])
    profile2.load(input_handles[1])

    diff.dynamic_smooth(profile1, profile2)

    profile1.save(output_handles[0])
    profile2.save(output_handles[1])
#smooth

def pair_diff(input_handles, distance_function="default", pairwise="diff-prod",
        pairwise_func="", do_smooth=False, summary="min", summary_func="",
        threshold=0, do_scale=False, down=False, do_positive=False,
        do_balance=False, precision=3):
    """
    Calculate the difference between two k-mer profiles.

    :arg input_handles: Open readable handles to a pair of k-mer profiles.
    :type input_handles: list(stream)
    :arg distance_function: Use a specific distance function.
    :type distance_function: str
    :arg pairwise: Name of the pairwise distance function.
    :type pairwise: str
    :arg pairwise_func: Custom pairwise distance function.
    :type pairwise_func: str
    :arg do_smooth: Enable smoothing.
    :type do_smooth: bool
    :arg summary: Name of the summary function.
    :type summary: str
    :arg summary_func: Custom summary function.
    :type summary_func: str
    :arg threshold: Threshold for the summary function.
    :type threshold: int
    :arg do_scale: Scale the profiles.
    :type do_scale: bool
    :arg down: Scale down.
    :type down: bool
    :arg do_positive: Only use positive values.
    :type do_positive: bool
    :arg do_balance: Balance the profiles.
    :type do_balance: bool
    :arg precision: Number of digits in the output.
    :type precision: int
    """
    summary_function = metrics.summary[summary]
    if summary_func:
        summary_function = eval("lambda " + summary_func)

    pairwise_function = metrics.pairwise[pairwise]
    if pairwise_func:
        pairwise_function = eval("lambda " + pairwise_func)

    diff = kdifflib.kMerDiff(do_balance=do_balance, do_positive=do_positive,
        do_smooth=do_smooth, summary=summary_function, threshold=threshold,
        do_scale=do_scale, down=down,
        distance_function=metrics.vector_distance[distance_function],
        pairwise=pairwise_function)

    profile1 = klib.kMer()
    profile2 = klib.kMer()

    profile1.load(input_handles[0])
    profile2.load(input_handles[1])

    if profile1.length != profile2.length:
        raise ValueError(length_error)

    print ("%%.%if" % precision) % diff.distance(profile1, profile2)
#pair_diff

def matrix_diff(input_handles, output_handle, distance_function="default",
        pairwise="diff-prod", pairwise_func="", do_smooth=False, summary="min",
        summary_func="", threshold=0, do_scale=False, down=False,
        do_positive=False, do_balance=False, precision=3):
    """
    Make a distance matrix any number of k-mer profiles.

    :arg input_handles: Open readable handles to a list of k-mer profiles.
    :type input_handles: list(stream)
    :arg output_handle: Open writeable handle to a distance matrix.
    :type output_handle: stream
    :arg euclidean: Use the Euclidean distance.
    :type euclidean: bool
    :arg pairwise: Name of the pairwise distance function.
    :type pairwise: str
    :arg pairwise_func: Custom pairwise distance function.
    :type pairwise_func: str
    :arg do_smooth: Enable smoothing.
    :type do_smooth: bool
    :arg summary: Name of the summary function.
    :type summary: str
    :arg summary_func: Custom summary function.
    :type summary_func: str
    :arg threshold: Threshold for the summary function.
    :type threshold: int
    :arg do_scale: Scale the profiles.
    :type do_scale: bool
    :arg down: Scale down.
    :type down: bool
    :arg do_positive: Only use positive values.
    :type do_positive: bool
    :arg do_balance: Balance the profiles.
    :type do_balance: bool
    :arg precision: Number of digits in the output.
    :type precision: int
    """
    if len(input_handles) < 2:
        raise ValueError("You must give at least two input files.")

    summary_function = metrics.summary[summary]
    if summary_func:
        summary_function = eval("lambda " + summary_func)

    pairwise_function = metrics.pairwise[pairwise]
    if pairwise_func:
        pairwise_function = eval("lambda " + pairwise_func)

    diff = kdifflib.kMerDiff(do_balance=do_balance, do_positive=do_positive,
        do_smooth=do_smooth, summary=summary_function, threshold=threshold,
        do_scale=do_scale, down=down,
        distance_function=metrics.vector_distance[distance_function],
        pairwise=pairwise_function)

    counts = []
    for i in input_handles:
        counts.append(klib.kMer())
        counts[-1].load(i)
        if counts[0].length != counts[-1].length:
            raise ValueError(length_error)
    #for

    kdifflib.distance_matrix(counts, output_handle, precision, diff)
#matrix_diff

def main():
    """
    Main entry point.
    """
    output_parser = argparse.ArgumentParser(add_help=False)
    output_parser.add_argument("output_handle", metavar="OUTPUT",
        type=ProtectedFileType('w'), help="output file")

    input_parser = argparse.ArgumentParser(add_help=False)
    input_parser.add_argument("input_handle", metavar="INPUT",
        type=argparse.FileType('r'), help="input file")

    pair_out_parser = argparse.ArgumentParser(add_help=False)
    pair_out_parser.add_argument("output_handles", metavar="OUTPUT", nargs=2,
        type=ProtectedFileType('w'), help="pair of output files")

    pair_in_parser = argparse.ArgumentParser(add_help=False)
    pair_in_parser.add_argument("input_handles", metavar="INPUT", nargs=2,
        type=argparse.FileType('r'), help="pair of input files")

    list_in_parser = argparse.ArgumentParser(add_help=False)
    list_in_parser.add_argument("input_handles", metavar="INPUT", nargs='+',
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
    diff_parser.add_argument("-b", dest="do_balance", default=False,
        action="store_true", help="balance the profiles (default=%(default)s)")
    diff_parser.add_argument("-p", dest="do_positive", default=False,
        action="store_true", help="use only positive values "
        "(default=%(default)s)")
    diff_parser.add_argument("-S", dest="do_scale", default=False,
        action="store_true", help="scale the profiles (default=%(default)s)")
    diff_parser.add_argument("-m", dest="do_smooth", default=False,
        action="store_true", help="smooth the profiles (default=%(default)s)")
    diff_parser.add_argument("-D", dest="distance_function", type=str,
        default="default", choices=metrics.vector_distance,
        help='choose distance function (%(type)s default="%(default)s")')
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
        output_parser], description=doc_split(index))
    parser_index.add_argument("size", metavar="SIZE", type=int,
        help="k-mer size (%(type)s)")
    parser_index.set_defaults(func=index)

    parser_merge = subparsers.add_parser("merge", parents=[pair_in_parser,
        output_parser], description=doc_split(merge))
    parser_merge.add_argument("-m", dest="merger", type=str,
        default="sum", choices=metrics.mergers,
        help='merge function (%(type)s default="%(default)s")')
    parser_merge.add_argument("--merge-function", dest="merge_func",
        type=str, default="", help="custom merge function "
        '(%(type)s default="%(default)s")')
    parser_merge.set_defaults(func=merge)

    parser_balance = subparsers.add_parser("balance", parents=[input_parser,
        output_parser], description=doc_split(balance))
    parser_balance.set_defaults(func=balance)

    parser_showbalance = subparsers.add_parser("showbalance",
        parents=[input_parser, precision_parser],
        description=doc_split(get_balance))
    parser_showbalance.set_defaults(func=get_balance)

    parser_stats = subparsers.add_parser("stats", 
        parents=[input_parser, precision_parser],
        description=doc_split(get_stats))
    parser_stats.set_defaults(func=get_stats)

    parser_distr = subparsers.add_parser("distr", 
        parents=[input_parser, output_parser],
        description=doc_split(distribution))
    parser_distr.set_defaults(func=distribution)

    parser_info = subparsers.add_parser("info", 
        parents=[input_parser],
        description=doc_split(info))
    parser_info.set_defaults(func=info)

    parser_getcount = subparsers.add_parser("getcount", 
        parents=[input_parser],
        description=doc_split(get_count))
    parser_getcount.add_argument("word", metavar="WORD", type=str,
        help="the word in question (%(type)s)")
    parser_getcount.set_defaults(func=get_count)

    parser_positive = subparsers.add_parser("positive",
        parents=[pair_in_parser, pair_out_parser],
        description=doc_split(positive))
    parser_positive.set_defaults(func=positive)

    parser_scale = subparsers.add_parser("scale", parents=[pair_in_parser,
        pair_out_parser, scale_parser], description=doc_split(scale))
    parser_scale.set_defaults(func=scale)

    parser_shrink = subparsers.add_parser("shrink", parents=[input_parser,
        output_parser], description=doc_split(shrink))
    parser_shrink.add_argument("factor", metavar="FACTOR", type=int,
        help="shrinking factor")
    parser_shrink.set_defaults(func=shrink)

    parser_shrink = subparsers.add_parser("shuffle", parents=[input_parser,
        output_parser], description=doc_split(shuffle))
    parser_shrink.set_defaults(func=shuffle)

    parser_smooth = subparsers.add_parser("smooth", parents=[pair_in_parser,
        pair_out_parser, smooth_parser], description=doc_split(smooth))
    parser_smooth.set_defaults(func=smooth)

    parser_diff = subparsers.add_parser("diff", parents=[diff_parser,
        pair_in_parser], description=doc_split(pair_diff))
    parser_diff.set_defaults(func=pair_diff)

    parser_matrix = subparsers.add_parser("matrix", parents=[diff_parser,
        list_in_parser, output_parser], description=doc_split(matrix_diff))
    parser_matrix.set_defaults(func=matrix_diff)

    try:
        arguments = parser.parse_args()
    except IOError, error:
        parser.error(error)

    try:
        arguments.func(**{k: v for k, v in vars(arguments).items()
            if k not in ("func", "subcommand")})
    except ValueError, error:
        parser.error(error)
#main

if __name__ == "__main__":
    main()
