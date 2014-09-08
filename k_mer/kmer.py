"""
Toolbox for k-mer profiles.
"""


from __future__ import division

import argparse
import os
import sys

from . import (ProtectedFileType, ProfileFileType, doc_split, usage, version,
               klib, kdifflib, metrics)

# Todo: Probably only used in user-defined custom functions, we should only
#   import it there.
from math import *


LENGTH_ERROR = 'k-mer lengths of the files differ'
NAMES_COUNT_ERROR = ('number of profile names does not match number of '
                     'profiles')
PAIRED_NAMES_COUNT_ERROR = ('number of left and right profile names do not '
                            'match')


def _name_from_handle(handle):
    """
    Try to get a name for `handle` from its filename, if there is one. Return
    `None` for non-filename handles such as string buffers or streams.
    """
    if not hasattr(handle, 'name') or handle.name.startswith('<'):
        return None
    return os.path.splitext(os.path.basename(handle.name))[0]


def convert(input_handles, output_handle, names=None):
    """
    Save k-mer profiles from files in the old plaintext format (used by kMer
    versions < 1.0.0) to a k-mer profile file in the current HDF5 format.

    :arg input_handles: Open readable k-mer profile file handles (old format).
    :type input_handles: list(file-like object)
    :arg output_handle: Open writeable k-mer profile file handle.
    :type output_handle: h5py.File
    :arg names: Optional list of names for the saved k-mer profiles (must have
      the same length as `input_handles`). If not provided, profiles are named
      according to the input filenames, or numbered consecutively from 1 if no
      filenames are available.
    :type names: list(str)
    """
    names = names or [_name_from_handle(h) for h in input_handles]

    if len(names) != len(input_handles):
        raise ValueError(NAMES_COUNT_ERROR)

    for input_handle, name in zip(input_handles, names):
        profile = klib.Profile.from_file_old_format(input_handle, name=name)
        profile.save(output_handle)


def index(input_handles, output_handle, size, names=None):
    """
    Make k-mer profiles from FASTA files.

    :arg input_handles: Open readable FASTA file handles.
    :type input_handles: list(file-like object)
    :arg output_handle: Open writeable k-mer profile file handle.
    :type output_handle: h5py.File
    :arg size: Size of k.
    :type size: int
    :arg names: Optional list of names for the created k-mer profiles (must
      have the same length as `input_handles`). If not provided, profiles are
      named according to the input filenames, or numbered consecutively from 1
      if no filenames are available.
    :type names: list(str)
    """
    names = names or [_name_from_handle(h) for h in input_handles]

    if len(names) != len(input_handles):
        raise ValueError(NAMES_COUNT_ERROR)

    for input_handle, name in zip(input_handles, names):
        profile = klib.Profile.from_fasta(input_handle, size, name=name)
        profile.save(output_handle)


def merge(input_handle_left, input_handle_right, output_handle,
          names_left=None, names_right=None, merger='sum', merge_func=None):
    """
    Merge k-mer profiles. If the files contain more than one profile, they are
    linked by name and merged pairwise. The resulting profile name is set to
    that of the original profiles if they match, or to their concatenation
    otherwise.

    :arg input_handle_left, input_handle_right: Open readable k-mer profile
      file handle.
    :type input_handle_left, input_handle_right: h5py.File
    :arg output_handle: Open writeable k-mer profile file handle.
    :type output_handle: h5py.File
    :arg names_left, names_right: Optional list of names of the k-mer profiles
      to consider. If not provided, all profiles in the file are considered.
    :type names_left, names_right: list(str)
    :arg merger: Merge function.
    :type merger: function
    :arg merge_func: Custom merge function.
    :type merge_func: str
    """
    names_left = names_left or sorted(input_handle_left['profiles'])
    names_right = names_right or sorted(input_handle_right['profiles'])

    if len(names_left) != len(names_right):
        raise ValueError(PAIRED_NAMES_COUNT_ERROR)

    if merge_func:
        # Todo: Use the more general approach from Wiggelen.
        # Todo: I guess this must be vectorized to an ufunc.
        merge_function = eval('lambda ' + merge_func)
    else:
        merge_function = metrics.mergers[merger]

    for name_left, name_right in zip(names_left, names_right):
        profile_left = klib.Profile.from_file(input_handle_left,
                                              name=name_left)
        profile_right = klib.Profile.from_file(input_handle_right,
                                               name=name_right)

        if profile_left.length != profile_right.length:
            raise ValueError(LENGTH_ERROR)

        profile_right.merge(profile_left, merge_function)

        if name_left == name_right:
            name = name_left
        else:
            name = name_left + '_' + name_right
        profile_right.save(output_handle, name=name)


def balance(input_handle, output_handle, names=None):
    """
    Balance k-mer profiles.

    :arg input_handle: Open readable k-mer profile file handle.
    :type input_handle: h5py.File
    :arg output_handle: Open writeable k-mer profile file handle.
    :type output_handle: h5py.File
    :arg names: Optional list of names of the k-mer profiles to consider. If
      not provided, all profiles in the file are considered in alphabetical
      order.
    :type names: list(str)
    """
    names = names or sorted(input_handle['profiles'])

    for name in names:
        profile = klib.Profile.from_file(input_handle, name=name)
        profile.balance()
        profile.save(output_handle)


def get_balance(input_handle, output_handle, precision=3, names=None):
    """
    Show the balance of k-mer profiles.

    :arg input_handle: Open readable k-mer profile file handle.
    :type input_handle: h5py.File
    :arg output_handle: Open writeable k-mer profile file handle.
    :type output_handle: file-like object
    :arg precision: Number of digits in the output.
    :type precision: int
    :arg names: Optional list of names of the k-mer profiles to consider. If
      not provided, all profiles in the file are considered in alphabetical
      order.
    :type names: list(str)
    """
    # Todo: Wouldn't it make more sense conceptually to create the reverse
    # complement profile and calculate the distance to that? Perhaps it's
    # even easier to vectorize in NumPy, so faster?
    names = names or sorted(input_handle['profiles'])

    for name in names:
        profile = klib.Profile.from_file(input_handle, name=name)

        forward, reverse = profile.split()
        output_handle.write(
            ('%%s %%.%if\n' % precision) %
            (name, metrics.multiset(forward, reverse,
                                    metrics.pairwise['diff-prod'])))


def get_stats(input_handle, output_handle, precision=3, names=None):
    """
    Show the mean and standard deviation of k-mer profiles.

    :arg input_handle: Open readable k-mer profile file handle.
    :type input_handle: h5py.File
    :arg output_handle: Open writeable file handle.
    :type output_handle: file-like object
    :arg precision: Number of digits in the output.
    :type precision: int
    :arg names: Optional list of names of the k-mer profiles to consider. If
      not provided, all profiles in the file are considered in alphabetical
      order.
    :type names: list(str)
    """
    names = names or sorted(input_handle['profiles'])

    for name in names:
        profile = klib.Profile.from_file(input_handle, name=name)

        output_handle.write(
            ('%%s %%.%if %%.%if\n' % (precision, precision)) %
            (name, profile.mean, profile.std))


def distribution(input_handle, output_handle, names=None):
    """
    Calculate the distribution of the values in k-mer profiles. Every output
    line has the name of the profile, the count and the number of k-mers with
    this count.

    :arg input_handle: Open readable k-mer profile file handle.
    :type input_handle: h5py.File
    :arg output_handle: Open writeable file handle.
    :type output_handle: file-like object
    :arg names: Optional list of names of the k-mer profiles to consider. If
      not provided, all profiles in the file are considered in alphabetical
      order.
    :type names: list(str)
    """
    names = names or sorted(input_handle['profiles'])

    for name in names:
        profile = klib.Profile.from_file(input_handle)

        output_handle.write(''.join('%s %i %i\n' % (name, v, c) for v, c in
                                    metrics.distribution(profile.counts)))


def info(input_handle, output_handle, names=None):
    """
    Print some information about k-mer profiles.

    :arg input_handle: Open readable k-mer profile file handle.
    :type input_handle: h5py.File
    :arg output_handle: Open writeable file handle.
    :type output_handle: file-like object
    :arg names: Optional list of names of the k-mer profiles to consider. If
      not provided, all profiles in the file are considered in alphabetical
      order.
    :type names: list(str)
    """
    names = names or sorted(input_handle['profiles'])

    output_handle.write('File format version: %s\n' %
                        input_handle.attrs['version'])
    output_handle.write('Produced by: %s\n' % input_handle.attrs['producer'])

    for name in names:
        # Todo: This loads the entire profile, which we actually do not need.
        #   Especially for files with many profiles, we might want to just
        #   get the statistics without loading the profile.
        profile = klib.Profile.from_file(input_handle, name=name)

        output_handle.write('\n' + '\n'.join([
            'Profile: %s' % profile.name,
            '- k-mer length: %i (%i k-mers)' % (profile.length,
                                                profile.number),
            '- Zero counts: %i' % (profile.number - profile.non_zero),
            '- Non-zero counts: %i' % profile.non_zero,
            '- Sum of counts: %i' % profile.total,
            '- Mean of counts: %.3f' % profile.mean,
            '- Median of counts: %i' % profile.median,
            '- Standard deviation of counts: %.3f' % profile.std]) + '\n')
        # Todo: Mean/median/std of non-zero counts.


def get_count(input_handle, output_handle, word, names=None):
    """
    Retrieve the counts in k-mer profiles for a particular word.

    :arg input_handle: Open readable handle to a k-mer profile file.
    :type input_handle: h5py.File
    :arg output_handle: Open writeable file handle.
    :type output_handle: file-like object
    :arg word: The query word.
    :type word: str
    :arg names: Optional list of names of the k-mer profiles to consider. If
      not provided, all profiles in the file are considered in alphabetical
      order.
    :type names: list(str)
    """
    names = names or sorted(input_handle['profiles'])

    for name in names:
        profile = klib.Profile.from_file(input_handle, name=name)
        if profile.length != len(word):
            raise ValueError('the length of the query does not match the '
                             'profile length')
        try:
            offset = profile.dna_to_binary(word)
        except KeyError:
            raise ValueError('the input is not a valid DNA sequence')
        output_handle.write('%s %i\n' % (name, profile.counts[offset]))


def positive(input_handle_left, input_handle_right, output_handle_left,
             output_handle_right, names_left=None, names_right=None):
    """
    Only keep counts that are positive in both k-mer profiles. If the files
    contain more than one profile, they are linked by name and processed
    pairwise.

    :arg input_handle_left, input_handle_right: Open readable k-mer profile
      file handle.
    :type input_handle_left, input_handle_right: h5py.File
    :arg output_handle_left, output_handle_right: Open writeable k-mer profile
      file handle.
    :type output_handle_left, output_handle_right: h5py.File
    :arg names_left, names_right: Optional list of names of the k-mer profiles
      to consider. If not provided, all profiles in the file are considered.
    :type names_left, names_right: list(str)
    """
    names_left = names_left or sorted(input_handle_left['profiles'])
    names_right = names_right or sorted(input_handle_right['profiles'])

    if len(names_left) != len(names_right):
        raise ValueError(PAIRED_NAMES_COUNT_ERROR)

    for name_left, name_right in zip(names_left, names_right):
        profile_left = klib.Profile.from_file(input_handle_left,
                                              name=name_left)
        profile_right = klib.Profile.from_file(input_handle_right,
                                               name=name_right)

        if profile_left.length != profile_right.length:
            raise ValueError(LENGTH_ERROR)

        profile_left.counts = metrics.positive(profile_left.counts,
                                               profile_right.counts)
        profile_right.counts = metrics.positive(profile_right.counts,
                                                profile_left.counts)

        profile_left.save(output_handle_left)
        profile_right.save(output_handle_right)


def scale(input_handle_left, input_handle_right, output_handle_left,
          output_handle_right, names_left=None, names_right=None, down=False):
    """
    Scale two profiles such that the total number of k-mers is equal. If the
    files contain more than one profile, they are linked by name and processed
    pairwise.

    :arg input_handle_left, input_handle_right: Open readable k-mer profile
      file handle.
    :type input_handle_left, input_handle_right: h5py.File
    :arg output_handle_left, output_handle_right: Open writeable k-mer profile
      file handle.
    :type output_handle_left, output_handle_right: h5py.File
    :arg names_left, names_right: Optional list of names of the k-mer profiles
      to consider. If not provided, all profiles in the file are considered.
    :type names_left, names_right: list(str)
    :arg down: Scale down.
    :type down: bool
    """
    names_left = names_left or sorted(input_handle_left['profiles'])
    names_right = names_right or sorted(input_handle_right['profiles'])

    if len(names_left) != len(names_right):
        raise ValueError(PAIRED_NAMES_COUNT_ERROR)

    for name_left, name_right in zip(names_left, names_right):
        profile_left = klib.Profile.from_file(input_handle_left,
                                              name=name_left)
        profile_right = klib.Profile.from_file(input_handle_right,
                                               name=name_right)

        if profile_left.length != profile_right.length:
            raise ValueError(LENGTH_ERROR)

        scale_left, scale_right = metrics.get_scale(profile_left.counts,
                                                    profile_right.counts)
        if down:
            scale_left, scale_right = metrics.scale_down(scale_left,
                                                         scale_right)
        profile_left.counts = profile_left.counts * scale_left
        profile_right.counts = profile_right.counts * scale_right

        profile_left.save(output_handle_left)
        profile_right.save(output_handle_right)


def shrink(input_handle, output_handle, factor, names=None):
    """
    Shrink k-mer profiles, effectively reducing k.

    :arg input_handle: Open readable k-mer profile file handle.
    :type input_handle: h5py.File
    :arg output_handle: Open writeable k-mer profile file handle.
    :type output_handle: h5py.File
    :arg factor: Scaling factor.
    :type factor: int
    :arg names: Optional list of names of the k-mer profiles to consider. If
      not provided, all profiles in the file are considered in alphabetical
      order.
    :type names: list(str)
    """
    names = names or sorted(input_handle['profiles'])

    for name in names:
        profile = klib.Profile.from_file(input_handle, name=name)
        profile.shrink(factor)
        profile.save(output_handle)


def shuffle(input_handle, output_handle, names=None):
    """
    Randomise k-mer profiles.

    :arg input_handle: Open readable k-mer profile file handle.
    :type input_handle: h5py.File
    :arg output_handle: Open writeable k-mer profile file handle.
    :type output_handle: h5py.File
    :arg names: Optional list of names of the k-mer profiles to consider. If
      not provided, all profiles in the file are considered in alphabetical
      order.
    :type names: list(str)
    """
    names = names or sorted(input_handle['profiles'])

    for name in names:
        profile = klib.Profile.from_file(input_handle, name=name)
        profile.shuffle()
        profile.save(output_handle)


def smooth(input_handle_left, input_handle_right, output_handle_left,
           output_handle_right, names_left=None, names_right=None,
           summary='min', summary_func=None, threshold=0):
    """
    Smooth two profiles by collapsing sub-profiles. If the files contain more
    than one profile, they are linked by name and processed pairwise.

    :arg input_handle_left, input_handle_right: Open readable k-mer profile
      file handle.
    :type input_handle_left, input_handle_right: h5py.File
    :arg output_handle_left, output_handle_right: Open writeable k-mer profile
      file handle.
    :type output_handle_left, output_handle_right: h5py.File
    :arg names_left, names_right: Optional list of names of the k-mer profiles
      to consider. If not provided, all profiles in the file are considered.
    :type names_left, names_right: list(str)
    :arg summary: Name of the summary function.
    :type summary: str
    :arg summary_func: Custom summary function.
    :type summary_func: str
    :arg threshold: Threshold for the summary function.
    :type threshold: int
    """
    names_left = names_left or sorted(input_handle_left['profiles'])
    names_right = names_right or sorted(input_handle_right['profiles'])

    if len(names_left) != len(names_right):
        raise ValueError(PAIRED_NAMES_COUNT_ERROR)

    if summary_func:
        # Todo: Use the more general approach from Wiggelen.
        # Todo: I guess this must be vectorized to an ufunc.
        smooth_function = eval('lambda ' + summary_func)
    else:
        smooth_function = metrics.summary[summary]

    diff = kdifflib.kMerDiff(summary=smooth_function, threshold=threshold)

    for name_left, name_right in zip(names_left, names_right):
        profile_left = klib.Profile.from_file(input_handle_left,
                                              name=name_left)
        profile_right = klib.Profile.from_file(input_handle_right,
                                               name=name_right)

        if profile_left.length != profile_right.length:
            raise ValueError(LENGTH_ERROR)

        diff.dynamic_smooth(profile_left, profile_right)

        profile_left.save(output_handle_left)
        profile_right.save(output_handle_right)


def pair_diff(input_handle_left, input_handle_right, output_handle,
              names_left=None, names_right=None, distance_function='default',
              pairwise='diff-prod', pairwise_func='', do_smooth=False,
              summary='min', summary_func='', threshold=0, do_scale=False,
              down=False, do_positive=False, do_balance=False, precision=3):
    """
    Calculate the difference between two k-mer profiles. If the files contain
    more than one profile, they are linked by name and processed pairwise.

    :arg input_handle_left, input_handle_right: Open readable k-mer profile
      file handle.
    :type input_handle_left, input_handle_right: h5py.File
    :arg output_handle: Open writeable file handle.
    :type output_handle: file-like object
    :arg names_left, names_right: Optional list of names of the k-mer profiles
      to consider. If not provided, all profiles in the file are considered.
    :type names_left, names_right: list(str)
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
        raise ValueError(PAIRED_NAMES_COUNT_ERROR)

    if summary_func:
        # Todo: Use the more general approach from Wiggelen.
        # Todo: I guess this must be vectorized to an ufunc.
        summary_function = eval('lambda ' + summary_func)
    else:
        summary_function = metrics.summary[summary]

    if pairwise_func:
        # Todo: Use the more general approach from Wiggelen.
        # Todo: I guess this must be vectorized to an ufunc.
        pairwise_function = eval('lambda ' + pairwise_func)
    else:
        pairwise_function = metrics.pairwise[pairwise]

    diff = kdifflib.kMerDiff(
        do_balance=do_balance, do_positive=do_positive, do_smooth=do_smooth,
        summary=summary_function, threshold=threshold, do_scale=do_scale,
        down=down, pairwise=pairwise_function,
        distance_function=metrics.vector_distance[distance_function])

    for name_left, name_right in zip(names_left, names_right):
        profile_left = klib.Profile.from_file(input_handle_left,
                                              name=name_left)
        profile_right = klib.Profile.from_file(input_handle_right,
                                               name=name_right)

        if profile_left.length != profile_right.length:
            raise ValueError(LENGTH_ERROR)

        output_handle.write(
            ('%%s %%s %%.%if\n' % precision) %
            (name_left, name_right,
             diff.distance(profile_left, profile_right)))


def matrix_diff(input_handle, output_handle, names=None,
                distance_function='default', pairwise='diff-prod',
                pairwise_func=None, do_smooth=False, summary='min',
                summary_func=None, threshold=0, do_scale=False, down=False,
                do_positive=False, do_balance=False, precision=3):
    """
    Make a distance matrix between any number of k-mer profiles.

    :arg input_handle: Open readable k-mer profile file handle.
    :type input_handle: h5py.File
    :arg output_handle: Open writeable file handle.
    :type output_handle: file-like object
    :arg names: Optional list of names of the k-mer profiles to consider. If
      not provided, all profiles in the file are considered in alphabetical
      order.
    :type names: list(str)
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
    names = names or sorted(input_handle['profiles'])

    if len(names) < 2:
        raise ValueError('you must give at least two k-mer profiles')

    if summary_func:
        # Todo: Use the more general approach from Wiggelen.
        # Todo: I guess this must be vectorized to an ufunc.
        summary_function = eval('lambda ' + summary_func)
    else:
        summary_function = metrics.summary[summary]

    if pairwise_func:
        # Todo: Use the more general approach from Wiggelen.
        # Todo: I guess this must be vectorized to an ufunc.
        pairwise_function = eval('lambda ' + pairwise_func)
    else:
        pairwise_function = metrics.pairwise[pairwise]

    diff = kdifflib.kMerDiff(
        do_balance=do_balance, do_positive=do_positive, do_smooth=do_smooth,
        summary=summary_function, threshold=threshold, do_scale=do_scale,
        down=down, pairwise=pairwise_function,
        distance_function=metrics.vector_distance[distance_function])

    # Todo: This first creates all profiles in memory, which may be
    #   problematic for large k. We could provide an option to keep only the
    #   two current profiles in memory. We may combine this with an option to
    #   paralellize this function. The downside is that profiles have to be
    #   read from file on each use.
    counts = []
    for name in names:
        counts.append(klib.Profile.from_file(input_handle, name=name))
        if counts[0].length != counts[-1].length:
            raise ValueError(LENGTH_ERROR)

    kdifflib.distance_matrix(counts, output_handle, precision, diff)


def main():
    """
    Command line interface.
    """
    multi_input_parser = argparse.ArgumentParser(add_help=False)
    multi_input_parser.add_argument(
        'input_handles', metavar='INPUT', type=argparse.FileType('r'),
        nargs='*', default=[sys.stdin], help='input file (default: stdin)')

    input_profile_parser = argparse.ArgumentParser(add_help=False)
    input_profile_parser.add_argument(
        'input_handle', metavar='INPUT', type=ProfileFileType('r'),
        help='input k-mer profile file',)
    input_profile_parser.add_argument(
        '-p', '--profiles', dest='names', metavar='NAME', nargs='+',
        help='names of the k-mer profiles to consider (default: all profiles '
        'in INPUT, in alphabetical order)')

    paired_input_profile_parser = argparse.ArgumentParser(add_help=False)
    paired_input_profile_parser.add_argument(
        'input_handle_left', metavar='INPUT_LEFT', type=ProfileFileType('r'),
        help='input k-mer profile file (left)')
    paired_input_profile_parser.add_argument(
        'input_handle_right', metavar='INPUT_RIGHT',
        type=ProfileFileType('r'), help='input k-mer profile file (right)')
    paired_input_profile_parser.add_argument(
        '-l', '--profiles-left', dest='names_left', metavar='NAME', nargs='+',
        help='names of the k-mer profiles to consider (left) (default: all '
        'profiles in INPUT_LEFT, in alphabetical order')
    paired_input_profile_parser.add_argument(
        '-r', '--profiles-right', dest='names_right', metavar='NAME',
        nargs='+', help='names of the k-mer profiles to consider (right) '
        '(default: all profiles in INPUT_RIGHT, in alphabetical order)')

    output_parser = argparse.ArgumentParser(add_help=False)
    output_parser.add_argument(
        'output_handle', metavar='OUTPUT', type=ProtectedFileType('w'),
        help='output file')

    output_profile_parser = argparse.ArgumentParser(add_help=False)
    output_profile_parser.add_argument(
        'output_handle', metavar='OUTPUT', type=ProfileFileType('w'),
        help='output k-mer profile file')

    paired_output_profile_parser = argparse.ArgumentParser(add_help=False)
    paired_output_profile_parser.add_argument(
        'output_handle_left', metavar='OUTPUT_LEFT',
        type=ProfileFileType('w'), help='output k-mer profile file (left)')
    paired_output_profile_parser.add_argument(
        'output_handle_right', metavar='OUTPUT_RIGHT',
        type=ProfileFileType('w'), help='output k-mer profile file (right)')

    scale_parser = argparse.ArgumentParser(add_help=False)
    scale_parser.add_argument(
        '-d', dest='down', action='store_true', help='scale down')

    smooth_parser = argparse.ArgumentParser(add_help=False)
    smooth_parser.add_argument(
        '-s', dest='summary', default='min', choices=metrics.summary,
        help='summary function for dynamic smoothing (default: %(default)s)')
    smooth_parser.add_argument(
        '--summary-function', metavar='STRING', dest='summary_func',
        help='custom summary function')
    smooth_parser.add_argument(
        '-t', dest='threshold', metavar='INT', type=int, default=0,
        help='threshold for the summary function (default: %(default)s)')

    precision_parser = argparse.ArgumentParser(add_help=False)
    precision_parser.add_argument(
        '-n', metavar='INT', dest='precision', type=int, default=3,
        help='precision in number of decimals (default: %(default)s)')

    diff_parser = argparse.ArgumentParser(
        add_help=False,
        parents=[scale_parser, smooth_parser, precision_parser])
    diff_parser.add_argument(
        '-b', '--balance', dest='do_balance', action='store_true',
        help='balance the profiles')
    diff_parser.add_argument(
        '--positive', dest='do_positive', action='store_true',
        help='use only positive values')
    diff_parser.add_argument(
        '-S', '--scale', dest='do_scale', action='store_true',
        help='scale the profiles')
    diff_parser.add_argument(
        '-m', '--smooth', dest='do_smooth', action='store_true',
        help='smooth the profiles')
    diff_parser.add_argument(
        '-D', dest='distance_function', default='default',
        choices=metrics.vector_distance,
        help='choose distance function (default: %(default)s)')
    diff_parser.add_argument(
        '-P', dest='pairwise', default='diff-prod', choices=metrics.pairwise,
        help='paiwise distance function for the multiset distance (default: '
        '%(default)s)')
    diff_parser.add_argument(
        '--pairwise-function', metavar='STRING', dest='pairwise_func',
        help='custom pairwise function')

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=usage[0], epilog=usage[1])
    parser.add_argument('-v', action="version", version=version(parser.prog))
    subparsers = parser.add_subparsers()

    parser_convert = subparsers.add_parser(
        'convert', parents=[output_profile_parser, multi_input_parser],
        description=doc_split(convert))
    parser_convert.add_argument(
        '-p', '--profiles', dest='names', metavar='NAME', nargs='+',
        help='names for the saved k-mer profiles, one per INPUT (default: '
        'profiles are named according to the input filenames, or numbered '
        'consecutively from 1 if no filenames are available)')
    parser_convert.set_defaults(func=convert)

    # Todo: Option to generate a profile per FASTA record instead of per FASTA
    #   file.
    parser_index = subparsers.add_parser(
        'index', parents=[output_profile_parser, multi_input_parser],
        description=doc_split(index))
    parser_index.add_argument(
        '-p', '--profiles', dest='names', metavar='NAME', nargs='+',
        help='names for the created k-mer profiles, one per INPUT (default: '
        'profiles are named according to the input filenames, or numbered '
        'consecutively from 1 if no filenames are available)')
    parser_index.add_argument(
        '-k', dest='size', metavar='SIZE', type=int, default=9,
        help='k-mer size (%(type)s default: %(default)s)')
    parser_index.set_defaults(func=index)

    parser_merge = subparsers.add_parser(
        'merge', parents=[paired_input_profile_parser, output_profile_parser],
        description=doc_split(merge))
    parser_merge.add_argument(
        '-m', dest='merger', type=str, default='sum',
        choices=metrics.mergers, help='merge function (default: %(default)s)')
    parser_merge.add_argument(
        '--merge-function', dest='merge_func', help='custom merge function')
    parser_merge.set_defaults(func=merge)

    parser_balance = subparsers.add_parser(
        'balance', parents=[input_profile_parser, output_profile_parser],
        description=doc_split(balance))
    parser_balance.set_defaults(func=balance)

    parser_showbalance = subparsers.add_parser(
        'showbalance', parents=[input_profile_parser, precision_parser],
        description=doc_split(get_balance))
    parser_showbalance.set_defaults(func=get_balance, output_handle=sys.stdout)

    parser_stats = subparsers.add_parser(
        'stats', parents=[input_profile_parser, precision_parser],
        description=doc_split(get_stats))
    parser_stats.set_defaults(func=get_stats, output_handle=sys.stdout)

    parser_distr = subparsers.add_parser(
        'distr', parents=[input_profile_parser, output_parser],
        description=doc_split(distribution))
    parser_distr.set_defaults(func=distribution)

    parser_info = subparsers.add_parser(
        'info', parents=[input_profile_parser],
        description=doc_split(info))
    parser_info.set_defaults(func=info, output_handle=sys.stdout)

    parser_getcount = subparsers.add_parser(
        'getcount', parents=[input_profile_parser],
        description=doc_split(get_count))
    parser_getcount.add_argument(
        'word', metavar='WORD', type=str, help='the word in question')
    parser_getcount.set_defaults(func=get_count, output_handle=sys.stdout)

    parser_positive = subparsers.add_parser(
        'positive', parents=[paired_input_profile_parser,
                             paired_output_profile_parser],
        description=doc_split(positive))
    parser_positive.set_defaults(func=positive)

    parser_scale = subparsers.add_parser(
        'scale', parents=[paired_input_profile_parser,
                          paired_output_profile_parser, scale_parser],
        description=doc_split(scale))
    parser_scale.set_defaults(func=scale)

    parser_shrink = subparsers.add_parser(
        'shrink', parents=[input_profile_parser, output_profile_parser],
        description=doc_split(shrink))
    parser_shrink.add_argument(
        '-f', '--factor', dest='factor', metavar='INT', type=int, default=1,
        help='shrinking factor (default: %(default)s)')
    parser_shrink.set_defaults(func=shrink)

    parser_shuffle = subparsers.add_parser(
        'shuffle', parents=[input_profile_parser, output_profile_parser],
        description=doc_split(shuffle))
    parser_shuffle.set_defaults(func=shuffle)

    parser_smooth = subparsers.add_parser(
        'smooth', parents=[paired_input_profile_parser,
                           paired_output_profile_parser, smooth_parser],
        description=doc_split(smooth))
    parser_smooth.set_defaults(func=smooth)

    parser_diff = subparsers.add_parser(
        'diff', parents=[paired_input_profile_parser, diff_parser],
        description=doc_split(pair_diff))
    parser_diff.set_defaults(func=pair_diff, output_handle=sys.stdout)

    parser_matrix = subparsers.add_parser(
        'matrix', parents=[input_profile_parser, output_parser, diff_parser],
        description=doc_split(matrix_diff))
    parser_matrix.set_defaults(func=matrix_diff)

    try:
        arguments = parser.parse_args()
    except IOError, error:
        parser.error(error)

    try:
        arguments.func(**dict((k, v) for k, v in vars(arguments).items()
                              if k not in ('func', 'subcommand')))
    except ValueError, error:
        parser.error(error)
