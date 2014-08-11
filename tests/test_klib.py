"""
Tests for the `k_mer.klib` module.
"""


import collections
import itertools
import os
import shutil
import tempfile

from k_mer import klib

try:
    # Python 2.7 and up.
    from collections import Counter
except ImportError:
    from _counter import Counter


# Some 8-mers.
LENGTH_8 = ['GTACATGA',
            'TAAACTAA',
            'TATCTTTA',
            'TACTATGT']

# Some 8-mers with N bases.
LENGTH_8_N = ['GNACATGA',
              'TAAACTNA',
              'TATNTTTA',
              'TACTATGN']

# Some 60-mers.
LENGTH_60 = ['GTACATGATAGGTCCACAGCTCTGAGCAAGGCAGACGTCCATACTTAAAACCCAGACTGC',
             'TAAACTAAAAGAAAGAATTTTTTTAATGGTAGACTACCTAAAATTATGTCTCTTAGTCCT',
             'TATCTTTACCTATATATTTGACTAAGATTTAGTATTACTACTACCTAAAATTATGTCTCT',
             'TACTATGTCTTGAAGGACAGCACCTGACCTCCCCCTGCAAGGTGTCATCCCCAAGCTGGT']

# Some 60-mers with N bases.
LENGTH_60_N = ['GTANATGATAGGTCCACAGCTCTGAGCAAGGCAGACGTCCATACTTAAAACCCAGACTGC',
               'TAAACTAAAAGAAAGAATTTTTTTAATGGTAGACTACCTAAAATTATGTCTCTTAGNCCT',
               'TATCTTTACCTATATATTTGACTAAGATTTAGTNTTACTACTACCTAAAATTATGTCTCT',
               'TACTATGTCTTGAAGGACAGCCCTGACCTCCCCCTGCAAGGTGTCATCCCCAAGCTGGTN']


def _fasta(sequences):
    """
    A serialization of `sequences` in FASTA format.
    """
    names = ('>sequence_%d' % i for i in itertools.count(1))
    return '\n'.join('\n'.join(entry) for entry in
                     itertools.izip(names, sequences)) + '\n'


def _count_index(sequence):
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


def _count(profile, sequence):
    """
    The count of `sequence` in a kMer profile.
    """
    return profile.count[_count_index(sequence)]


def _counts(sequence, k):
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


class TestKlib():
    def setup(self):
        self.temp_dir = tempfile.mkdtemp(prefix='kmer-tests-')

    def teardown(self):
        shutil.rmtree(self.temp_dir)

    def _touch(self, content=None):
        """
        Create a file and optionally write `contents` to it. Filename is
        returned.
        """
        os_handle, filename = tempfile.mkstemp(dir=self.temp_dir)
        os.close(os_handle)
        if content is not None:
            with open(filename, 'w') as f:
                f.write(content)
        return filename

    def _test_profile(self, profile, counts, k):
        """
        Validate the given kMer profile, using `counts` as a reference.
        """
        assert profile.length == k
        assert profile.total == sum(counts.values())
        assert profile.non_zero == len(counts)
        for s in counts:
            assert _count(profile, s) == counts[s]

    def _test_analyse(self, sequences, k):
        """
        Test `k_mer.klib.analyse` with `sequences` and `k=k`.
        """
        counts = _counts(sequences, k)
        profile = klib.kMer()
        with open(self._touch(_fasta(sequences))) as fasta_handle:
            profile.analyse(fasta_handle, k)
        self._test_profile(profile, counts, k)

    def test_analyse_single_k1(self):
        self._test_analyse(LENGTH_60[:1], 1)

    def test_analyse_single_k4(self):
        self._test_analyse(LENGTH_60[:1], 4)

    def test_analyse_single_strlen(self):
        self._test_analyse(LENGTH_8[:1], 8)

    def test_analyse_single_almost_strlen(self):
        self._test_analyse(LENGTH_8[:1], 7)

    def test_analyse_multi_k1(self):
        self._test_analyse(LENGTH_60, 1)

    def test_analyse_multi_k4(self):
        self._test_analyse(LENGTH_60, 4)

    def test_analyse_multi_strlen(self):
        self._test_analyse(LENGTH_8, 8)

    def test_analyse_multi_almost_strlen(self):
        self._test_analyse(LENGTH_8, 7)

    def test_analyse_multi_n_k1(self):
        self._test_analyse(LENGTH_60_N, 1)

    def test_analyse_multi_n_k4(self):
        self._test_analyse(LENGTH_60_N, 4)

    def test_analyse_multi_n_strlen(self):
        self._test_analyse(LENGTH_8_N, 8)

    def test_analyse_multi_n_almost_strlen(self):
        self._test_analyse(LENGTH_8_N, 7)

    def test_profile_save(self):
        profile = klib.kMer()
        with open(self._touch(_fasta(LENGTH_60))) as fasta_handle:
            profile.analyse(fasta_handle, 4)

        filename = self._touch()
        with open(filename, 'w') as profile_handle:
            profile.save(profile_handle)

        counts = _counts(LENGTH_60, 4)
        with open(filename) as profile_handle:
            assert int(next(profile_handle)) == 4
            assert int(next(profile_handle)) == sum(counts.values())
            assert int(next(profile_handle)) == len(counts)
            profile_counts = [int(l) for l in profile_handle]
            for s in counts:
                assert profile_counts[_count_index(s)] == counts[s]

    def test_profile_load(self):
        counts = _counts(LENGTH_60, 4)
        content = '4\n%d\n%d\n' % (sum(counts.values()), len(counts))
        content += '\n'.join(str(counts[''.join(s)]) for s in
                             itertools.product('ACGT', repeat=4)) + '\n'
        filename = self._touch(content)

        profile = klib.kMer()
        with open(filename) as profile_handle:
            profile.load(profile_handle)
        assert profile.length == 4
        assert profile.total == sum(counts.values())
        assert profile.non_zero == len(counts)
        for s in counts:
            assert _count(profile, s) == counts[s]

    def test_profile_load_save(self):
        counts = _counts(LENGTH_60, 4)
        content = '4\n%d\n%d\n' % (sum(counts.values()), len(counts))
        content += '\n'.join(str(counts[''.join(s)]) for s in
                             itertools.product('ACGT', repeat=4)) + '\n'
        filename_load = self._touch(content)
        filename_save = self._touch()

        profile = klib.kMer()
        with open(filename_load) as profile_handle:
            profile.load(profile_handle)

        with open(filename_save, 'w') as profile_handle:
            profile.save(profile_handle)

        with open(filename_save) as f:
            assert f.read() == content

    def test_profile_save_load(self):
        profile_save = klib.kMer()
        with open(self._touch(_fasta(LENGTH_60))) as fasta_handle:
            profile_save.analyse(fasta_handle, 4)

        filename = self._touch()
        with open(filename, 'w') as profile_handle:
            profile_save.save(profile_handle)

        profile_load = klib.kMer()
        with open(filename) as profile_handle:
            profile_load.load(profile_handle)

        assert profile_save.length == profile_load.length
        assert profile_save.total == profile_load.total
        assert profile_save.non_zero == profile_load.non_zero
        assert profile_save.count == profile_load.count
