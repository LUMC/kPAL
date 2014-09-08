"""
Tests for the `k_mer.klib` module.
"""


from __future__ import division

import itertools

from Bio import Seq
import numpy as np
import pytest

from k_mer import klib

import utils


class TestKlib(utils.TestEnvironment):
    def test_profile(self):
        counts = utils.counts(utils.SEQUENCES, 8)
        profile = klib.Profile(utils.as_array(counts, 8))
        utils.test_profile(profile, counts, 8)

    def _test_from_fasta(self, sequences, k):
        counts = utils.counts(sequences, k)
        with open(self.fasta(sequences)) as fasta_handle:
            profile = klib.Profile.from_fasta(fasta_handle, k)
        utils.test_profile(profile, counts, k)

    def test_from_fasta_single_k1(self):
        self._test_from_fasta(utils.SEQUENCES[:1], 1)

    def test_from_fasta_single_k4(self):
        self._test_from_fasta(utils.SEQUENCES[:1], 4)

    def test_from_fasta_single_strlen(self):
        self._test_from_fasta(utils.LENGTH_8[:1], 8)

    def test_from_fasta_single_almost_strlen(self):
        self._test_from_fasta(utils.LENGTH_8[:1], 7)

    def test_from_fasta_multi_k1(self):
        self._test_from_fasta(utils.SEQUENCES, 1)

    def test_from_fasta_multi_k4(self):
        self._test_from_fasta(utils.SEQUENCES, 4)

    def test_from_fasta_multi_strlen(self):
        self._test_from_fasta(utils.LENGTH_8, 8)

    def test_from_fasta_multi_almost_strlen(self):
        self._test_from_fasta(utils.LENGTH_8, 7)

    def test_from_fasta_multi_n_k1(self):
        self._test_from_fasta(utils.SEQUENCES_WITH_N, 1)

    def test_from_fasta_multi_n_k4(self):
        self._test_from_fasta(utils.SEQUENCES_WITH_N, 4)

    def test_from_fasta_multi_n_strlen(self):
        self._test_from_fasta(utils.LENGTH_8_WITH_N, 8)

    def test_from_fasta_multi_n_almost_strlen(self):
        self._test_from_fasta(utils.LENGTH_8_WITH_N, 7)

    def test_profile_save(self):
        counts = utils.counts(utils.SEQUENCES, 4)
        profile = klib.Profile(utils.as_array(counts, 4))

        filename = self.empty()
        with utils.open_profile(filename, 'w') as profile_handle:
            profile.save(profile_handle)

        utils.test_profile_file(filename, counts, 4)

    def test_profile_from_file_old_format(self):
        counts = utils.counts(utils.SEQUENCES, 4)
        with open(self.profile_old_format(counts, 4)) as handle:
            profile = klib.Profile.from_file_old_format(handle)

        utils.test_profile(profile, counts, 4)

    def test_profile_from_file(self):
        counts = utils.counts(utils.SEQUENCES, 4)
        with utils.open_profile(self.profile(counts, 4), 'r') as profile_handle:
            profile = klib.Profile.from_file(profile_handle)

        utils.test_profile(profile, counts, 4)

    def test_profile_from_file_save(self):
        counts = utils.counts(utils.SEQUENCES, 4)
        with utils.open_profile(self.profile(counts, 4), 'r') as profile_handle:
            profile = klib.Profile.from_file(profile_handle)

        filename = self.empty()
        with utils.open_profile(filename, 'w') as profile_handle:
            profile.save(profile_handle)

        utils.test_profile_file(filename, counts, 4)

    def test_profile_save_from_file(self):
        counts = utils.counts(utils.SEQUENCES, 4)
        profile = klib.Profile(utils.as_array(counts, 4))

        filename = self.empty()
        with utils.open_profile(filename, 'w') as profile_handle:
            profile.save(profile_handle)

        with utils.open_profile(filename, 'r') as profile_handle:
            profile = klib.Profile.from_file(profile_handle)

        utils.test_profile(profile, counts, 4)

    def test_profile_merge(self):
        counts_left = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_right = utils.counts(utils.SEQUENCES_RIGHT, 8)

        profile_left = klib.Profile(utils.as_array(counts_left, 8))
        profile_right = klib.Profile(utils.as_array(counts_right, 8))

        profile_left.merge(profile_right)
        utils.test_profile(profile_left, counts_left + counts_right, 8)

    def test_profile_balance(self):
        counts = utils.counts(utils.SEQUENCES, 8)
        profile = klib.Profile(utils.as_array(counts, 8))
        profile.balance()

        counts.update(dict((str(Seq.reverse_complement(s)), c)
                           for s, c in counts.items()))
        utils.test_profile(profile, counts, 8)

    def test_profile_balance_palindrome(self):
        counts = utils.counts(['AATT'], 4)
        profile = klib.Profile(utils.as_array(counts, 4))
        profile.balance()

        counts.update(dict((str(Seq.reverse_complement(s)), c)
                           for s, c in counts.items()))
        utils.test_profile(profile, counts, 4)

    def _test_profile_split(self, sequences, length):
        counts = utils.counts(sequences, length)
        profile = klib.Profile(utils.as_array(counts, length))
        left, right = profile.split()

        assert len(left) == len(right)
        assert sum(left) + sum(right) == sum(counts.values()) * 2

        indices_left = {}
        indices_right = {}
        indices_palindrome = {}

        for s, c in counts.items():
            r = str(Seq.reverse_complement(s))
            if s < r:
                indices_left[utils.count_index(s)] = c * 2
            elif s > r:
                indices_right[utils.count_index(r)] = counts[s] * 2
            else:
                indices_palindrome[utils.count_index(s)] = c

        assert ([c for c in left if c > 0] ==
                [c for i, c in sorted(indices_left.items() +
                                      indices_palindrome.items())])
        assert ([c for c in right if c > 0] ==
                [c for i, c in sorted(indices_right.items() +
                                      indices_palindrome.items())])

    def test_profile_split(self):
        self._test_profile_split(utils.SEQUENCES, 8)

    def test_profile_split_short(self):
        self._test_profile_split(utils.SEQUENCES, 2)

    def test_profile_split_palindrome(self):
        self._test_profile_split(utils.SEQUENCES + ['ACCTAGGT'], 8)

    def test_profile_reverse_complement(self):
        counts = utils.counts(utils.SEQUENCES, 8)
        profile = klib.Profile(utils.as_array(counts, 8))

        for i in range(profile.length):
            assert (profile.binary_to_dna(profile.reverse_complement(i)) ==
                    str(Seq.reverse_complement(profile.binary_to_dna(i))))

    def test_profile_reverse_complement_palindrome(self):
        counts = utils.counts(['ACCTAGGT'], 8)
        profile = klib.Profile(utils.as_array(counts, 8))

        for i in range(profile.length):
            assert (profile.binary_to_dna(profile.reverse_complement(i)) ==
                    str(Seq.reverse_complement(profile.binary_to_dna(i))))

    def test_profile_shrink(self):
        counts = utils.counts(utils.SEQUENCES, 8)
        profile = klib.Profile(utils.as_array(counts, 8))
        profile.shrink(1)

        counts = utils.Counter(dict((t, sum(counts[u] for u in counts
                                            if u.startswith(t)))
                                    for t in set(s[:-1] for s in counts)))
        utils.test_profile(profile, counts, 7)

    def test_profile_shrink_two(self):
        counts = utils.counts(utils.SEQUENCES, 8)
        profile = klib.Profile(utils.as_array(counts, 8))
        profile.shrink(2)

        counts = utils.Counter(dict((t, sum(counts[u] for u in counts
                                            if u.startswith(t)))
                                    for t in set(s[:-2] for s in counts)))
        utils.test_profile(profile, counts, 6)

    def test_profile_shrink_three(self):
        counts = utils.counts(utils.SEQUENCES, 8)
        profile = klib.Profile(utils.as_array(counts, 8))
        profile.shrink(3)

        counts = utils.Counter(dict((t, sum(counts[u] for u in counts
                                            if u.startswith(t)))
                                    for t in set(s[:-3] for s in counts)))
        utils.test_profile(profile, counts, 5)

    def test_profile_shrink_max(self):
        counts = utils.counts(utils.SEQUENCES, 4)
        profile = klib.Profile(utils.as_array(counts, 4))
        profile.shrink(3)

        counts = utils.Counter(dict((t, sum(counts[u] for u in counts
                                            if u.startswith(t)))
                                    for t in set(s[:-3] for s in counts)))
        utils.test_profile(profile, counts, 1)

    def test_profile_shrink_invalid(self):
        counts = utils.counts(utils.SEQUENCES, 4)
        profile = klib.Profile(utils.as_array(counts, 4))
        with pytest.raises(ValueError):
            profile.shrink(4)

    def test_profile_shuffle(self):
        counts = utils.counts(utils.SEQUENCES, 2)
        profile = klib.Profile(utils.as_array(counts, 2))

        np.random.seed(100)
        profile.shuffle()

        counts = dict(zip([''.join(s) for s in itertools.product('ACGT', repeat=2)],
                          [13,  7,  6, 18, 12,  1, 13, 17, 16, 12, 23, 27, 24, 17, 18, 12]))
        utils.test_profile(profile, counts, 2)

    def test_profile_dna_to_binary(self):
        counts = utils.counts(utils.SEQUENCES, 4)
        profile = klib.Profile(utils.as_array(counts, 4))

        for i, s in enumerate(itertools.product('ACGT', repeat=4)):
            assert i == profile.dna_to_binary(''.join(s))

    def test_profile_binary_to_dna(self):
        counts = utils.counts(utils.SEQUENCES, 4)
        profile = klib.Profile(utils.as_array(counts, 4))

        for i, s in enumerate(itertools.product('ACGT', repeat=4)):
            assert ''.join(s) == profile.binary_to_dna(i)

    def test_profile_ratios_matrix(self):
        counts = utils.counts(utils.SEQUENCES, 4)
        profile = klib.Profile(utils.as_array(counts, 4))

        ratios = profile._ratios_matrix()
        total = sum(counts.values())

        for left in itertools.product('ACGT', repeat=4):
            for right in itertools.product('ACGT', repeat=4):
                left = ''.join(left)
                right = ''.join(right)
                ratio = ratios[utils.count_index(left)][utils.count_index(right)]
                try:
                    assert ratio == counts[left] / counts[right] / total
                except ZeroDivisionError:
                    assert ratio == -1.0

    def test_profile_freq_diff_matrix(self):
        counts = utils.counts(utils.SEQUENCES, 4)
        profile = klib.Profile(utils.as_array(counts, 4))
        freq_diffs = profile._freq_diff_matrix()

        total = sum(counts.values())

        for left in itertools.product('ACGT', repeat=4):
            for right in itertools.product('ACGT', repeat=4):
                left = ''.join(left)
                right = ''.join(right)
                freq_diff = freq_diffs[utils.count_index(left)][utils.count_index(right)]
                if counts[right] > 0:
                    assert freq_diff == abs(counts[left] - counts[right]) / total
                else:
                    assert freq_diff == 0

    def test_profile_print_counts(self, capsys):
        counts = utils.counts(utils.SEQUENCES, 4)
        profile = klib.Profile(utils.as_array(counts, 4))
        profile.print_counts()

        out, err = capsys.readouterr()
        assert out == ''.join('%s %d\n' % (''.join(s), counts[''.join(s)])
                              for s in itertools.product('ACGT', repeat=4))
