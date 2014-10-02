"""
Toolbox for k-mer profiles.

.. moduleauthor:: Leiden University Medical Center <humgen@lumc.nl>
.. moduleauthor:: Jeroen F.J. Laros <j.f.j.laros@lumc.nl>
.. moduleauthor:: Martijn Vermaat <m.vermaat.hg@lumc.nl>

.. Licensed under the MIT license, see the LICENSE file.
"""


from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import str, zip

import argparse
import importlib
import os
import re
import sys

import numpy as np

from . import (USAGE, FileType, ProfileFileType, doc_split, version, klib,
               kdistlib, metrics)


LENGTH_ERROR = 'k-mer lengths of the files differ'
NAMES_COUNT_ERROR = ('number of profile names does not match number of '
                     'profiles')
PAIRED_NAMES_COUNT_ERROR = ('number of left and right profile names do not '
                            'match')


# Importable definition, e.g. `package.module.merge_function`.
# http://docs.python.org/2/reference/lexical_analysis.html#identifiers
_PYTHON_IMPORTABLE = '{0}(\.{0})+$'.format('[_a-zA-Z][_a-zA-Z0-9]*')


def _name_from_handle(handle):
    """
    Try to get a name for `handle` from its filename, if there is one. Return
    `None` for non-filename handles such as string buffers or streams.
    """
    if not hasattr(handle, 'name') or handle.name.startswith('<'):
        return None
    return os.path.splitext(os.path.basename(str(handle.name)))[0]


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


def count(input_handles, output_handle, size, names=None):
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
          names_left=None, names_right=None, merger='sum', custom_merger=None):
    """
    Merge k-mer profiles. If the files contain more than one profile, they are
    linked by name and merged pairwise. The resulting profile name is set to
    that of the original profiles if they match, or to their concatenation
    otherwise.

    :arg h5py.File input_handle_left, input_handle_right: Open readable k-mer
      profile file handle.
    :arg h5py.File output_handle: Open writeable k-mer profile file handle.
    :arg list(str) names_left, names_right: Optional list of names of the
      k-mer profiles to consider. If not provided, all profiles in the file
      are considered.
    :arg function merger: Merge function.
    :arg str custom_merger: Custom merge function.
    """
    names_left = names_left or sorted(input_handle_left['profiles'])
    names_right = names_right or sorted(input_handle_right['profiles'])

    if len(names_left) != len(names_right):
        raise ValueError(PAIRED_NAMES_COUNT_ERROR)

    if custom_merger:
        if re.match(_PYTHON_IMPORTABLE, custom_merger):
            # Importable definition, e.g. `package.module.merge_function`.
            module, name = custom_merger.rsplit('.', 1)
            merge_function = getattr(importlib.import_module(module), name)
        else:
            # Expression over `left` and `right`, e.g. `np.add(left, right)`.
            # The `numpy` package is available as `np`.
            merge_function = eval('lambda left, right: ' + custom_merger,
                                  {'np': np})
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

    :arg h5py.File input_handle: Open readable k-mer profile file handle.
    :arg h5py.File output_handle: Open writeable k-mer profile file handle.
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

    :arg h5py.File input_handle: Open readable k-mer profile file handle.
    :arg output_handle: Open writeable k-mer profile file handle.
    :type output_handle: file-like object
    :arg int precision: Number of digits in the output.
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
        balance = metrics.multiset(forward, reverse,
                                   metrics.pairwise['prod'])
        print(name, '{{0:.{0}f}}'.format(precision).format(balance),
              file=output_handle)


def get_stats(input_handle, output_handle, precision=3, names=None):
    """
    Show the mean and standard deviation of k-mer profiles.

    :arg h5py.File input_handle: Open readable k-mer profile file handle.
    :arg output_handle: Open writeable file handle.
    :type output_handle: file-like object
    :arg int precision: Number of digits in the output.
    :arg names: Optional list of names of the k-mer profiles to consider. If
      not provided, all profiles in the file are considered in alphabetical
      order.
    :type names: list(str)
    """
    names = names or sorted(input_handle['profiles'])

    for name in names:
        profile = klib.Profile.from_file(input_handle, name=name)

        print(name,
              '{{0:.{0}f}}'.format(precision).format(profile.mean),
              '{{0:.{0}f}}'.format(precision).format(profile.std),
              file=output_handle)


def distribution(input_handle, output_handle, names=None):
    """
    Calculate the distribution of the values in k-mer profiles. Every output
    line has the name of the profile, the count and the number of k-mers with
    this count.

    :arg h5py.File input_handle: Open readable k-mer profile file handle.
    :arg output_handle: Open writeable file handle.
    :type output_handle: file-like object
    :arg names: Optional list of names of the k-mer profiles to consider. If
      not provided, all profiles in the file are considered in alphabetical
      order.
    :type names: list(str)
    """
    names = names or sorted(input_handle['profiles'])

    for name in names:
        profile = klib.Profile.from_file(input_handle, name=name)

        print('\n'.join('{0} {1} {2}'.format(name, v, c)
                        for v, c in metrics.distribution(profile.counts)),
              file=output_handle)


def info(input_handle, output_handle, names=None):
    """
    Print some information about k-mer profiles.

    :arg h5py.File input_handle: Open readable k-mer profile file handle.
    :arg output_handle: Open writeable file handle.
    :type output_handle: file-like object
    :arg names: Optional list of names of the k-mer profiles to consider. If
      not provided, all profiles in the file are considered in alphabetical
      order.
    :type names: list(str)
    """
    names = names or sorted(input_handle['profiles'])

    print('File format version:', input_handle.attrs['version'],
          file=output_handle)
    print('Produced by:', input_handle.attrs['producer'], file=output_handle)

    for name in names:
        # Todo: This loads the entire profile, which we actually do not need.
        #   Especially for files with many profiles, we might want to just
        #   get the statistics without loading the profile.
        profile = klib.Profile.from_file(input_handle, name=name)

        print('', file=output_handle)
        print('Profile:', profile.name, file=output_handle)
        print('- k-mer length:', str(profile.length),
              '({0} k-mers)'.format(profile.number), file=output_handle)
        print('- Zero counts:', str(profile.number - profile.non_zero),
              file=output_handle)
        print('- Non-zero counts:', str(profile.non_zero), file=output_handle)
        print('- Sum of counts:', str(profile.total), file=output_handle)
        print('- Mean of counts:', '{0:.3f}'.format(profile.mean),
              file=output_handle)
        print('- Median of counts:', '{0:.3f}'.format(profile.median),
              file=output_handle)
        print('- Standard deviation of counts:',
              '{0:.3f}'.format(profile.std), file=output_handle)


def get_count(input_handle, output_handle, word, names=None):
    """
    Retrieve the counts in k-mer profiles for a particular word.

    :arg h5py.File input_handle: Open readable handle to a k-mer profile file.
    :arg output_handle: Open writeable file handle.
    :type output_handle: file-like object
    :arg str word: The query word.
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
        print(name, str(profile.counts[offset]), file=output_handle)


def positive(input_handle_left, input_handle_right, output_handle_left,
             output_handle_right, names_left=None, names_right=None):
    """
    Only keep counts that are positive in both k-mer profiles. If the files
    contain more than one profile, they are linked by name and processed
    pairwise.

    :arg h5py.File input_handle_left, input_handle_right: Open readable k-mer
      profile file handle.
    :arg h5py.File output_handle_left, output_handle_right: Open writeable
      k-mer profile file handle.
    :arg list(str) names_left, names_right: Optional list of names of the
      k-mer profiles consider. If not provided, all profiles in the file are
      considered.
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

    :arg h5py.File input_handle_left, input_handle_right: Open readable k-mer
      profile handle.
    :arg h5py.File output_handle_left, output_handle_right: Open writeable
      k-mer profile handle.
    :arg list(str) names_left, names_right: Optional list of names of the
      k-mer profiles to consider. If not provided, all profiles in the file
      are considered.
    :arg bool down: Scale down.
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

    :arg h5py.File input_handle: Open readable k-mer profile file handle.
    :arg h5py.File output_handle: Open writeable k-mer profile file handle.
    :arg int factor: Scaling factor.
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

    :arg h5py.File input_handle: Open readable k-mer profile file handle.
    :arg h5py.File output_handle: Open writeable k-mer profile file handle.
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
           summary='min', custom_summary=None, threshold=0):
    """
    Smooth two profiles by collapsing sub-profiles. If the files contain more
    than one profile, they are linked by name and processed pairwise.

    :arg h5py.File input_handle_left, input_handle_right: Open readable k-mer
      profile handle.
    :arg h5py.File output_handle_left, output_handle_right: Open writeable
      k-mer profile file handle.
    :arg list(str) names_left, names_right: Optional list of names of the
      k-mer profiles to consider. If not provided, all profiles in the file
      are considered.
    :arg str summary: Name of the summary function.
    :arg str custom_summary: Custom summary function.
    :arg int threshold: Threshold for the summary function.
    """
    names_left = names_left or sorted(input_handle_left['profiles'])
    names_right = names_right or sorted(input_handle_right['profiles'])

    if len(names_left) != len(names_right):
        raise ValueError(PAIRED_NAMES_COUNT_ERROR)

    if custom_summary:
        if re.match(_PYTHON_IMPORTABLE, custom_summary):
            # Importable definition, e.g. `package.module.summary_function`.
            module, name = custom_summary.rsplit('.', 1)
            summary_function = getattr(importlib.import_module(module), name)
        else:
            # Expression over `values`, e.g. `np.max(values)`. The `numpy`
            # package is available as `np`.
            summary_function = eval('lambda values: ' + custom_summary,
                                    {'np': np})
    else:
        summary_function = metrics.summary[summary]

    dist = kdistlib.ProfileDistance(summary=summary_function,
                                    threshold=threshold)

    for name_left, name_right in zip(names_left, names_right):
        profile_left = klib.Profile.from_file(input_handle_left,
                                              name=name_left)
        profile_right = klib.Profile.from_file(input_handle_right,
                                               name=name_right)

        if profile_left.length != profile_right.length:
            raise ValueError(LENGTH_ERROR)

        dist.dynamic_smooth(profile_left, profile_right)

        profile_left.save(output_handle_left)
        profile_right.save(output_handle_right)


def distance(input_handle_left, input_handle_right, output_handle,
             names_left=None, names_right=None, distance_function='default',
             pairwise='prod', custom_pairwise=None, do_smooth=False,
             summary='min', custom_summary=None, threshold=0, do_scale=False,
             down=False, do_positive=False, do_balance=False, precision=3):
    """
    Calculate the distance between two k-mer profiles. If the files contain
    more than one profile, they are linked by name and processed pairwise.

    :arg h5py.File input_handle_left, input_handle_right: Open readable k-mer
      profile handle.
    :arg output_handle: Open writeable file handle.
    :type output_handle: file-like object
    :arg list(str) names_left, names_right: Optional list of names of the
      k-mer profiles to consider. If not provided, all profiles in the file
      are considered.
    :arg str distance_function: Use a specific distance function.
    :arg str pairwise: Name of the pairwise distance function.
    :arg str custom_pairwise: Custom pairwise distance function.
    :arg bool do_smooth: Enable smoothing.
    :arg str summary: Name of the summary function.
    :arg str custom_summary: Custom summary function.
    :arg int threshold: Threshold for the summary function.
    :arg bool do_scale: Scale the profiles.
    :arg bool down: Scale down.
    :arg bool do_positive: Only use positive values.
    :arg bool do_balance: Balance the profiles.
    :arg int precision: Number of digits in the output.
    """
    names_left = names_left or sorted(input_handle_left['profiles'])
    names_right = names_right or sorted(input_handle_right['profiles'])

    if len(names_left) != len(names_right):
        raise ValueError(PAIRED_NAMES_COUNT_ERROR)

    if custom_summary:
        if re.match(_PYTHON_IMPORTABLE, custom_summary):
            # Importable definition, e.g. `package.module.summary_function`.
            module, name = custom_summary.rsplit('.', 1)
            summary_function = getattr(importlib.import_module(module), name)
        else:
            # Expression over `values`, e.g. `np.max(values)`. The `numpy`
            # package is available as `np`.
            summary_function = eval('lambda values: ' + custom_summary,
                                    {'np': np})
    else:
        summary_function = metrics.summary[summary]

    if custom_pairwise:
        if re.match(_PYTHON_IMPORTABLE, custom_pairwise):
            # Importable definition, e.g. `package.module.pairwise_function`.
            module, name = custom_pairwise.rsplit('.', 1)
            pairwise_function = getattr(importlib.import_module(module), name)
        else:
            # Expression over `left` and `right`, e.g. `np.add(left, right)`.
            # The `numpy` package is available as `np`.
            pairwise_function = eval('lambda left, right: ' + custom_pairwise,
                                     {'np': np})
    else:
        pairwise_function = metrics.pairwise[pairwise]

    dist = kdistlib.ProfileDistance(
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

        print(name_left, name_right,
              '{{0:.{0}f}}'.format(precision).format(dist.distance(
                  profile_left, profile_right)),
              file=output_handle)


def distance_matrix(input_handle, output_handle, names=None,
                    distance_function='default', pairwise='prod',
                    custom_pairwise=None, do_smooth=False, summary='min',
                    custom_summary=None, threshold=0, do_scale=False,
                    down=False, do_positive=False, do_balance=False,
                    precision=3):
    """
    Make a distance matrix between any number of k-mer profiles.

    :arg h5py.File input_handle: Open readable k-mer profile file handle.
    :arg output_handle: Open writeable file handle.
    :type output_handle: file-like object
    :arg names: Optional list of names of the k-mer profiles to consider. If
      not provided, all profiles in the file are considered in alphabetical
      order.
    :type names: list(str)
    :arg bool euclidean: Use the Euclidean distance.
    :arg str pairwise: Name of the pairwise distance function.
    :arg str custom_pairwise: Custom pairwise distance function.
    :arg bool do_smooth: Enable smoothing.
    :arg str summary: Name of the summary function.
    :arg str custom_summary: Custom summary function.
    :arg int threshold: Threshold for the summary function.
    :arg bool do_scale: Scale the profiles.
    :arg bool down: Scale down.
    :arg bool do_positive: Only use positive values.
    :arg bool do_balance: Balance the profiles.
    :arg int precision: Number of digits in the output.
    """
    names = names or sorted(input_handle['profiles'])

    if len(names) < 2:
        raise ValueError('you must give at least two k-mer profiles')

    if custom_summary:
        if re.match(_PYTHON_IMPORTABLE, custom_summary):
            # Importable definition, e.g. `package.module.summary_function`.
            module, name = custom_summary.rsplit('.', 1)
            summary_function = getattr(importlib.import_module(module), name)
        else:
            # Expression over `values`, e.g. `np.max(values)`. The `numpy`
            # package is available as `np`.
            summary_function = eval('lambda values: ' + custom_summary,
                                    {'np': np})
    else:
        summary_function = metrics.summary[summary]

    if custom_pairwise:
        if re.match(_PYTHON_IMPORTABLE, custom_pairwise):
            # Importable definition, e.g. `package.module.pairwise_function`.
            module, name = custom_pairwise.rsplit('.', 1)
            pairwise_function = getattr(importlib.import_module(module), name)
        else:
            # Expression over `left` and `right`, e.g. `np.add(left, right)`.
            # The `numpy` package is available as `np`.
            pairwise_function = eval('lambda left, right: ' + custom_pairwise,
                                     {'np': np})
    else:
        pairwise_function = metrics.pairwise[pairwise]

    dist = kdistlib.ProfileDistance(
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

    kdistlib.distance_matrix(counts, output_handle, precision, dist)


def main():
    """
    Command line interface.
    """
    multi_input_parser = argparse.ArgumentParser(add_help=False)
    multi_input_parser.add_argument(
        'input_handles', metavar='INPUT', type=FileType('r'), nargs='*',
        default=[sys.stdin], help='input file (default: stdin)')

    input_profile_parser = argparse.ArgumentParser(add_help=False)
    input_profile_parser.add_argument(
        'input_handle', metavar='INPUT', type=ProfileFileType('r'),
        help='input k-mer profile file',)
    input_profile_parser.add_argument(
        '-p', '--profiles', dest='names', metavar='NAME', type=str, nargs='+',
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
        '-l', '--profiles-left', dest='names_left', metavar='NAME', type=str,
        nargs='+', help='names of the k-mer profiles to consider (left) '
        '(default: all profiles in INPUT_LEFT, in alphabetical order')
    paired_input_profile_parser.add_argument(
        '-r', '--profiles-right', dest='names_right', metavar='NAME',
        type=str, nargs='+', help='names of the k-mer profiles to consider '
        '(right) (default: all profiles in INPUT_RIGHT, in alphabetical '
        'order)')

    output_parser = argparse.ArgumentParser(add_help=False)
    output_parser.add_argument(
        'output_handle', metavar='OUTPUT', type=FileType('w'),
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
        '-s', dest='summary', type=str, default='min',
        choices=metrics.summary, help='summary function for dynamic '
        'smoothing (default: %(default)s)')
    smooth_parser.add_argument(
        '-M', '--custom-summary', metavar='STRING', type=str,
        dest='custom_summary',
        help='custom Python summary function, specified either by an '
        'expression over the NumPy ndarray "values" (e.g., '
        '"np.max(values)"), or an importable name (e.g., '
        '"package.module.summary") that can be called with an ndarray as '
        'argument')
    smooth_parser.add_argument(
        '-t', dest='threshold', metavar='INT', type=int, default=0,
        help='threshold for the summary function (default: %(default)s)')

    precision_parser = argparse.ArgumentParser(add_help=False)
    precision_parser.add_argument(
        '-n', metavar='INT', dest='precision', type=int, default=3,
        help='precision in number of decimals (default: %(default)s)')

    dist_parser = argparse.ArgumentParser(
        add_help=False,
        parents=[scale_parser, smooth_parser, precision_parser])
    dist_parser.add_argument(
        '-b', '--balance', dest='do_balance', action='store_true',
        help='balance the profiles')
    dist_parser.add_argument(
        '--positive', dest='do_positive', action='store_true',
        help='use only positive values')
    dist_parser.add_argument(
        '-S', '--scale', dest='do_scale', action='store_true',
        help='scale the profiles')
    dist_parser.add_argument(
        '-m', '--smooth', dest='do_smooth', action='store_true',
        help='smooth the profiles')
    dist_parser.add_argument(
        '-D', dest='distance_function', type=str, default='default',
        choices=metrics.vector_distance,
        help='choose distance function (default: %(default)s)')
    dist_parser.add_argument(
        '-P', dest='pairwise', type=str, default='prod',
        choices=metrics.pairwise, help='paiwise distance function for the '
        'multiset distance (default: %(default)s)')
    # Todo: Note in the documentation that the custom pairwise must be
    #   vectorized with a tip for using np.vectorize if it is not.
    dist_parser.add_argument(
        '-f', '--pairwise-function', metavar='STRING', dest='custom_pairwise',
        type=str,
        help='custom Python pairwise function, specified either by an '
        'expression over the two NumPy ndarrays "left" and "right" (e.g., '
        '"abs(left - right) / (left + right + 1)"), or an importable name '
        '(e.g., "package.module.pairwise") that can be called with two '
        'ndarrays as arguments')

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=USAGE[0], epilog=USAGE[1])
    parser.add_argument('-v', action='version', version=version(parser.prog))
    subparsers = parser.add_subparsers()

    parser_convert = subparsers.add_parser(
        'convert', parents=[multi_input_parser, output_profile_parser],
        description=doc_split(convert))
    parser_convert.add_argument(
        '-p', '--profiles', dest='names', metavar='NAME', type=str, nargs='+',
        help='names for the saved k-mer profiles, one per INPUT (default: '
        'profiles are named according to the input filenames, or numbered '
        'consecutively from 1 if no filenames are available)')
    parser_convert.set_defaults(func=convert)

    # Todo: Option to generate a profile per FASTA record instead of per FASTA
    #   file.
    parser_count = subparsers.add_parser(
        'count', parents=[multi_input_parser, output_profile_parser],
        description=doc_split(count))
    parser_count.add_argument(
        '-p', '--profiles', dest='names', metavar='NAME', type=str, nargs='+',
        help='names for the created k-mer profiles, one per INPUT (default: '
        'profiles are named according to the input filenames, or numbered '
        'consecutively from 1 if no filenames are available)')
    parser_count.add_argument(
        '-k', dest='size', metavar='SIZE', type=int, default=9,
        help='k-mer size (%(type)s default: %(default)s)')
    parser_count.set_defaults(func=count)

    parser_merge = subparsers.add_parser(
        'merge', parents=[paired_input_profile_parser, output_profile_parser],
        description=doc_split(merge))
    parser_merge.add_argument(
        '-m', dest='merger', type=str, default='sum',
        choices=metrics.mergers, help='merge function (default: %(default)s)')
    # Todo: Note in the documentation that the custom merger must be
    #   vectorized with a tip for using np.vectorize if it is not.
    parser_merge.add_argument(
        '-c', '--custom-merger', dest='custom_merger', metavar='STRING',
        type=str,
        help='custom Python merge function, specified either by an '
        'expression over the two NumPy ndarrays "left" and "right" (e.g., '
        '"np.add(left, right)"), or an importable name (e.g., '
        '"package.module.merge") that can be called with two ndarrays as '
        'arguments')
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

    parser_distance = subparsers.add_parser(
        'distance', parents=[paired_input_profile_parser, dist_parser],
        description=doc_split(distance))
    parser_distance.set_defaults(func=distance, output_handle=sys.stdout)

    # Todo: I think we should just write to stdout.
    parser_matrix = subparsers.add_parser(
        'matrix', parents=[input_profile_parser, output_parser, dist_parser],
        description=doc_split(distance_matrix))
    parser_matrix.set_defaults(func=distance_matrix)

    try:
        arguments = parser.parse_args()
    except IOError as error:
        parser.error(error)

    try:
        arguments.func(**dict((k, v) for k, v in vars(arguments).items()
                              if k not in ('func', 'subcommand')))
    except ValueError as error:
        parser.error(error)
