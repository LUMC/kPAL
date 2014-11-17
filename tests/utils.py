"""
Utilities for kPAL unit tests.
"""


from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import range, str, zip
from future import standard_library

from io import open
import itertools
import os
import shutil
import tempfile

import h5py
import numpy as np

with standard_library.hooks():
    from collections import Counter


# Some 8-mers.
LENGTH_8 = ['GTACATGA',
            'TAAACTAA',
            'TATCTTTA',
            'TACTATGT']

# Some 8-mers with N bases.
LENGTH_8_WITH_N = ['GNACATGA',
                   'TAAACTNA',
                   'TATNTTTA',
                   'TACTATGN']

# Some 60-mers.
LENGTH_60 = ['GTACATGATAGGTCCACAGCTCTGAGCAAGGCAGACGTCCATACTTAAAACCCAGACTGC',
             'TAAACTAAAAGAAAGAATTTTTTTAATGGTAGACTACCTAAAATTATGTCTCTTAGTCCT',
             'TATCTTTACCTATATATTTGACTAAGATTTAGTATTACTACTACCTAAAATTATGTCTCT',
             'TACTATGTCTTGAAGGACAGCACCTGACCTCCCCCTGCAAGGTGTCATCCCCAAGCTGGT']

# Some 60-mers with N bases.
LENGTH_60_WITH_N = ['GTANATGATAGGTCCACAGCTCTGAGCAAGGCAGACGTCCATACTTAAAACCCAGACTGC',
                    'TAAACTAAAAGAAAGAATTTTTTTAATGGTAGACTACCTAAAATTATGTCTCTTAGNCCT',
                    'TATCTTTACCTATATATTTGACTAAGATTTAGTNTTACTACTACCTAAAATTATGTCTCT',
                    'TACTATGTCTTGAAGGACAGCCCTGACCTCCCCCTGCAAGGTGTCATCCCCAAGCTGGTN']

# Some more 60-mers.
LENGTH_60_MORE = ['TTACAATGATTAGGTCCACAGCTCTGAGCAACGCGCAGACGTCACATACTTCAAAACCCA',
                  'TAAAAACTATATAAGAAATCGAATTTTCTCTTAATGGTAGCAGCTACCGTAAAATCTATG',
                  'TACGCCTATATATCTTTGACTAAGCATTTATGTATTACATACTAACCAAAATTACTGTCT',
                  'TACTAAGTTTTGAAGGACAGCACATCACCTGCGCAATATCGGTGTCACCCCATAGCTCCT']

# Some alternative names for these k-mers.
SEQUENCES = LENGTH_60
SEQUENCES_SHORT = LENGTH_8
SEQUENCES_LEFT, SEQUENCES_RIGHT = LENGTH_60, LENGTH_60_MORE
SEQUENCES_WITH_N = LENGTH_60_WITH_N
SEQUENCES_SHORT_WITH_N = LENGTH_8_WITH_N


def reverse_complement(sequence):
    """
    Reverse complement of a sequence represented as unicode string.
    """
    complement = dict(zip([ord(b) for b in u'ACGTacgt'], u'TGCAtgca'))
    return ''.join(reversed(sequence.translate(complement)))


def counts(sequence, k):
    """
    Simple and naive k-mer counting. Returns a dictionary of `k`-mers with
    their counts in `sequence` (implemented as `collections.Counter`).

    To be used as a reference.
    """
    if isinstance(sequence, str):
        if sequence.startswith('>'):
            sequences = [s.split('\n')[1] for s in sequence.split('>')[1:]]
        else:
            sequences = [sequence]
    else:
        sequences = sequence

    counts = Counter()
    for sequence in sequences:
        for i in range(0, len(sequence) - k + 1):
            kmer = sequence[i:i + k]
            if 'N' not in kmer and 'n' not in kmer:
                counts[kmer] += 1
    return counts


def as_array(counts, k):
    """
    Convert dictionary of k-mer counts to a list containing the count for
    every possible k-mer.
    """
    return np.array([counts[''.join(s)] for s in
                     itertools.product('ACGT', repeat=k)])


def count_index(sequence):
    """
    The index of `sequence` in a kMer profile.
    """
    nucleotide_to_binary = {
        'A': 0x00, 'a': 0x00,
        'C': 0x01, 'c': 0x01,
        'G': 0x02, 'g': 0x02,
        'T': 0x03, 't': 0x03
    }
    binary = 0x00
    for b in sequence:
        binary = ((binary << 2) | nucleotide_to_binary[b])
    return binary


def count(profile, sequence):
    """
    The count of `sequence` in a kMer profile.
    """
    return profile.counts[count_index(sequence)]


class open_profile(h5py.File):
    """
    Context manager for an open kMer profile file.
    """
    def __init__(self, filename, mode='r'):
        super(open_profile, self).__init__(filename, mode=mode)
        if 'w' in mode:
            self.attrs['format'] = 'kMer'
            self.attrs['version'] = '1.0.0'
            self.attrs['producer'] = 'kMer unit tests'
            self.create_group('profiles')


def test_profile(profile, counts, k):
    """
    Validate the given kMer profile, using `counts` as a reference.
    """
    assert profile.length == k
    assert profile.total == sum(counts.values())
    assert profile.non_zero == len(counts)
    assert np.array_equal(profile.counts, as_array(counts, k))


def test_profile_file(filename, counts, k, name=None):
    """
    Validate the kMer in the given filename, using `counts` as a reference.
    """
    with open_profile(filename) as f:
        name = name or list(f['profiles'].keys())[0]
        profile = f['/profiles/%s' % name]
        assert profile.attrs['length'] == k
        assert profile.attrs['total'] == sum(counts.values())
        assert profile.attrs['non_zero'] == len(counts)
        assert np.array_equal(profile[:], as_array(counts, k))


class TestEnvironment(object):
    """
    Test class providing an isolated test environment for each test.
    """
    def setup(self):
        self.temp_dir = tempfile.mkdtemp(prefix='kmer-tests-')

    def teardown(self):
        shutil.rmtree(self.temp_dir)

    def empty(self):
        """
        Create an empty file and return filename.
        """
        os_handle, filename = tempfile.mkstemp(dir=self.temp_dir)
        os.close(os_handle)
        return filename

    def fasta(self, sequences):
        """
        Create a file and write `sequences` to it in FASTA format. Filename is
        returned.
        """
        filename = self.empty()

        with open(filename, 'w') as f:
            names = ('>sequence_%d' % i for i in itertools.count(1))
            f.write('\n'.join('\n'.join(entry) for entry in
                              zip(names, sequences)) + '\n')

        return filename

    def profile_old_format(self, counts, k):
        """
        Create a file and write `counts` to it in the old plaintext format.
        Filename is returned.
        """
        content = '%d\n%d\n%d\n' % (k, sum(counts.values()), len(counts))
        content += '\n'.join(str(counts[''.join(s)]) for s in
                             itertools.product('ACGT', repeat=k)) + '\n'

        filename = self.empty()

        with open(filename, 'w') as f:
            f.write(content)

        return filename

    def profile(self, counts, k, name='1'):
        """
        Create a file and write `counts` to it in HDF5 format. Filename is
        returned.
        """
        filename = self.empty()

        with open_profile(filename, 'w') as f:
            profile = f.create_dataset(
                'profiles/%s' % name, dtype='int64', compression='gzip',
                data=[counts[''.join(s)] for s in itertools.product('ACGT', repeat=k)])
            profile.attrs['length'] = k
            profile.attrs['total'] = sum(counts.values())
            profile.attrs['non_zero'] = len(counts)

        return filename

    def multi_profile(self, k, counts_list, names=None):
        """
        Create a file and write multiple `counts` to it in HDF5 format.
        Filename is returned.
        """
        names = names or [str(i + 1) for i in range(len(counts_list))]

        filename = self.empty()

        with open_profile(filename, 'w') as f:
            for counts, name in zip(counts_list, names):
                profile = f.create_dataset(
                    'profiles/%s' % name, dtype='int64', compression='gzip',
                    data=[counts[''.join(s)]
                          for s in itertools.product('ACGT', repeat=k)])
                profile.attrs['length'] = k
                profile.attrs['total'] = sum(counts.values())
                profile.attrs['non_zero'] = len(counts)

        return filename
