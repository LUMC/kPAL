#!/usr/bin/python

"""
Toolbox for k-mer profiles.
"""

from __future__ import division

import argparse
import itertools
import os
import sys

from math import *

from . import (ProtectedFileType, ProfileFileType, doc_split, usage, version,
               klib, kdifflib, metrics)

length_error = "k-mer lengths of the files differ."
names_count_error = "number of profile names does not match number of profiles"
pair_names_count_error = "number of left and right profile names do not match"

def name_from_handle(handle):
    if not hasattr(handle, 'name') or handle.name.startswith('<'):
        return None
    return os.path.splitext(os.path.basename(handle.name))[0]

def index(input_handles, output_handle, size, names=None):
    """
    Make a k-mer profile from a FASTA file.

    :arg input_handle: Open readable handle to a FASTA file.
    :type input_handle: stream
    :arg output_handle: Open writeable handle to a k-mer profile.
    :type output_handle: stream
    :arg size: Size of k.
    :type size: int
    """
    names = names or [name_from_handle(h) for h in input_handles]

    if len(names) != len(input_handles):
        raise ValueError(names_count_error)

    for input_handle, name in zip(input_handles, names):
        profile = klib.kMer.from_fasta(input_handle, size, name=name)
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

    profile1 = klib.kMer.from_file(input_handles[0])
    profile2 = klib.kMer.from_file(input_handles[1])

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
    profile = klib.kMer.from_file(input_handle)
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
    # Todo: Wouldn't it make more sense conceptually to create the reverse
    # complement profile and calculate the distance to that?
    profile = klib.kMer.from_file(input_handle)

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
    profile = klib.kMer.from_file(input_handle)

    print ("%%.%if %%.%if" % (precision, precision)) % \
        metrics.stats(profile.counts)
#get_stats

def distribution(input_handle, output_handle):
    """
    Calculate the distribution of the values in a k-mer profile.

    :arg input_handle: Open readable handle to a k-mer profile.
    :type input_handle: stream
    :arg output_handle: Open writeable handle to distribution file.
    :type output_handle: stream
    """
    profile = klib.kMer.from_file(input_handle)

    for i in metrics.distribution(profile.counts):
        output_handle.write("%i %i\n" % i)
#distribution

def info(input_handle, names=None):
    """
    Print some information about the k-mer profile.

    :arg input_handle: Open readable handle to a k-mer profile.
    :type input_handle: stream
    """
    names = names or sorted(input_handle['profiles'])

    print "File format version: %s" % input_handle.attrs['version']
    print "Produced by: %s" % input_handle.attrs['producer']

    for name in names:
        # Todo: This loads the entire profile, which we actually do not need.
        #   Especially for files with many profiles, we might want to just
        #   get the statistics without loading the profile.
        profile = klib.kMer.from_file(input_handle, name=name)
        print
        print "Profile: %s" % profile.name
        print "k-mer length: %i" % profile.length
        print "Total number of counts: %i" % profile.total
        print "Non-zero counts: %i" % profile.non_zero
#info

def get_count(input_handle, word, names=None):
    """
    Retrieve the count for a particular word.

    :arg input_handle: Open readable handle to a k-mer profile.
    :type input_handle: stream
    :arg word: The query word.
    :type word: str
    """
    names = names or sorted(input_handle['profiles'])

    for name in names:
        profile = klib.kMer.from_file(input_handle, name=name)
        if profile.length != len(word):
            raise ValueError("The length of the query does not match the profile "
                             "length.")
        try:
            offset = profile.dna_to_binary(word)
        except KeyError:
            raise ValueError("The input is not a valid DNA sequence.")
        print name, profile.counts[offset]
#get_count

def positive(input_handles, output_handles):
    """
    Only keep counts that are positive in both profiles.

    :arg input_handles: Open readable handles to a pair of k-mer profiles.
    :type input_handles: list(stream)
    :arg output_handles: Open writeable handles to a pair of k-mer profiles.
    :type output_handles: list(stream)
    """
    profile1 = klib.kMer.from_file(input_handles[0])
    profile2 = klib.kMer.from_file(input_handles[1])

    profile1.counts = metrics.positive(profile1.counts, profile2.counts)
    profile2.counts = metrics.positive(profile2.counts, profile1.counts)

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
    profile1 = klib.kMer.from_file(input_handles[0])
    profile2 = klib.kMer.from_file(input_handles[1])

    scale1, scale2 = metrics.get_scale(profile1.counts, profile2.counts)
    if down:
        scale1, scale2 = metrics.scale_down(scale1, scale2)
    profile1.counts = metrics.scale(profile1.counts, scale1)
    profile2.counts = metrics.scale(profile2.counts, scale2)

    profile1.save(output_handles[0])
    profile2.save(output_handles[1])
#scale

def shrink(input_handle, output_handle, factor, names=None):
    """
    Shrink a profile, effectively reducing k.

    :arg input_handle: Open readable handle to a k-mer profile.
    :type input_handle: stream
    :arg output_handle: Open writeable handle to a k-mer profile.
    :type output_handle: stream
    :arg factor: Scaling factor.
    :type factor: int
    """
    names = names or sorted(input_handle['profiles'])

    for name in names:
        profile = klib.kMer.from_file(input_handle, name=name)
        profile.shrink(factor)
        profile.save(output_handle)
#shrink

def shuffle(input_handle, output_handle, names=None):
    """
    Randomise a profile.

    :arg input_handle: Open readable handle to a k-mer profile.
    :type input_handle: stream
    :arg output_handle: Open writeable handle to a k-mer profile.
    :type output_handle: stream
    """
    names = names or sorted(input_handle['profiles'])

    for name in names:
        profile = klib.kMer.from_file(input_handle, name=name)
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
    profile1 = klib.kMer.from_file(input_handles[0])
    profile2 = klib.kMer.from_file(input_handles[1])

    diff.dynamic_smooth(profile1, profile2)

    profile1.save(output_handles[0])
    profile2.save(output_handles[1])
#smooth

def pair_diff(input_handle_left, input_handle_right, names_left=None, names_right=None,
        distance_function="default", pairwise="diff-prod", pairwise_func="",
        do_smooth=False, summary="min", summary_func="", threshold=0,
        do_scale=False, down=False, do_positive=False, do_balance=False,
        precision=3):
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
    names_left = names_left or sorted(input_handle_left['profiles'])
    names_right = names_right or sorted(input_handle_right['profiles'])

    if len(names_left) != len(names_right):
        raise ValueError(pair_names_count_error)

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

    for name_left, name_right in zip(names_left, names_right):
        profile_left = klib.kMer.from_file(input_handle_left, name=name_left)
        profile_right = klib.kMer.from_file(input_handle_right, name=name_right)

        if profile_left.length != profile_right.length:
            raise ValueError(length_error)

        print ("%%s %%s %%.%if" % precision) % (
            name_left, name_right, diff.distance(profile_left, profile_right))
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
        counts.append(klib.kMer.from_file(i))
        if counts[0].length != counts[-1].length:
            raise ValueError(length_error)
    #for

    kdifflib.distance_matrix(counts, output_handle, precision, diff)
#matrix_diff

def main():
    """
    Main entry point.
    """
    names_parser = argparse.ArgumentParser(add_help=False)
    names_parser.add_argument("-n", "--name", dest="names", metavar="NAME",
        nargs='+', help="profile name(s)")

    pair_names_parser = argparse.ArgumentParser(add_help=False)
    pair_names_parser.add_argument("-nl", "--name-left", dest="names_left",
        metavar="NAME", nargs='+', help="profile name(s) (left)")
    pair_names_parser.add_argument("-nr", "--name-right", dest="names_right",
        metavar="NAME", nargs='+', help="profile name(s) (right)")

    output_parser = argparse.ArgumentParser(add_help=False)
    output_parser.add_argument("output_handle", metavar="OUTPUT",
        type=ProtectedFileType('w'), help="output file")

    output_profile_parser = argparse.ArgumentParser(add_help=False)
    output_profile_parser.add_argument("output_handle", metavar="OUTPUT",
        type=ProfileFileType('w'), help="output file")

    input_parser = argparse.ArgumentParser(add_help=False)
    input_parser.add_argument("input_handle", metavar="INPUT",
        type=argparse.FileType('r'), help="input file")

    multi_input_parser = argparse.ArgumentParser(add_help=False)
    multi_input_parser.add_argument("input_handles", metavar="INPUT",
        type=argparse.FileType('r'), nargs='*', default=[sys.stdin],
        help="input file(s) (default=stdin)")

    input_profile_parser = argparse.ArgumentParser(add_help=False)
    input_profile_parser.add_argument("input_handle", metavar="INPUT",
        type=ProfileFileType('r'), help="input file")

    pair_input_profile_parser = argparse.ArgumentParser(add_help=False)
    pair_input_profile_parser.add_argument("input_handle_left", metavar="INPUT_LEFT",
        type=ProfileFileType('r'), help="input file (left)")
    pair_input_profile_parser.add_argument("input_handle_right", metavar="INPUT_RIGHT",
        type=ProfileFileType('r'), help="input file (right)")

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

    # Todo: Instead of creating a new profile file, append profiles to an
    # existing file.
    parser_index = subparsers.add_parser("index", parents=[output_profile_parser,
        multi_input_parser, names_parser], description=doc_split(index))
    parser_index.add_argument("-k", dest="size", metavar="SIZE", type=int,
        default=9, help='k-mer size (%(type)s default=%(default)s)')
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
        parents=[input_profile_parser, names_parser],
        description=doc_split(info))
    parser_info.set_defaults(func=info)

    parser_getcount = subparsers.add_parser("getcount",
        parents=[input_profile_parser, names_parser],
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

    parser_shrink = subparsers.add_parser("shrink",
        parents=[input_profile_parser, output_profile_parser, names_parser],
        description=doc_split(shrink))
    parser_shrink.add_argument("-f", "--factor", dest="factor",
        metavar="FACTOR", type=int, default=1,
        help="shrinking factor (%(type)s default=%(default)s)")
    parser_shrink.set_defaults(func=shrink)

    parser_shuffle = subparsers.add_parser("shuffle",
        parents=[input_profile_parser, output_profile_parser, names_parser],
        description=doc_split(shuffle))
    parser_shuffle.set_defaults(func=shuffle)

    parser_smooth = subparsers.add_parser("smooth", parents=[pair_in_parser,
        pair_out_parser, smooth_parser], description=doc_split(smooth))
    parser_smooth.set_defaults(func=smooth)

    parser_diff = subparsers.add_parser("diff", parents=[diff_parser,
        pair_input_profile_parser, pair_names_parser],
        description=doc_split(pair_diff))
    parser_diff.set_defaults(func=pair_diff)

    parser_matrix = subparsers.add_parser("matrix", parents=[diff_parser,
        list_in_parser, output_parser], description=doc_split(matrix_diff))
    parser_matrix.set_defaults(func=matrix_diff)

    try:
        arguments = parser.parse_args()
    except IOError, error:
        parser.error(error)

    try:
        arguments.func(**dict((k, v) for k, v in vars(arguments).items()
                              if k not in ("func", "subcommand")))
    except ValueError, error:
        parser.error(error)
#main

if __name__ == "__main__":
    main()
