"""
Tests for the `k_mer.klib` module.
"""


from Bio import Seq
import pytest

from k_mer import klib

import utils


class TestKlib(utils.TestEnvironment):
    def _test_from_fasta(self, sequences, k):
        """
        Test `k_mer.klib.analyse` with `sequences` and `k=k`.
        """
        counts = utils.counts(sequences, k)
        with open(self.fasta(sequences)) as fasta_handle:
            profile = klib.kMer.from_fasta(fasta_handle, k)
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
        profile = klib.kMer(utils.as_array(counts, 4))

        filename = self.empty()
        with utils.open_profile(filename, 'w') as profile_handle:
            profile.save(profile_handle)

        utils.test_profile_file(filename, counts, 4)

    def test_profile_from_file(self):
        counts = utils.counts(utils.SEQUENCES, 4)
        with utils.open_profile(self.profile(counts, 4), 'r') as profile_handle:
            profile = klib.kMer.from_file(profile_handle)

        utils.test_profile(profile, counts, 4)

    def test_profile_from_file_save(self):
        counts = utils.counts(utils.SEQUENCES, 4)
        with utils.open_profile(self.profile(counts, 4), 'r') as profile_handle:
            profile = klib.kMer.from_file(profile_handle)

        filename = self.empty()
        with utils.open_profile(filename, 'w') as profile_handle:
            profile.save(profile_handle)

        utils.test_profile_file(filename, counts, 4)

    def test_profile_save_from_file(self):
        counts = utils.counts(utils.SEQUENCES, 4)
        profile = klib.kMer(utils.as_array(counts, 4))

        filename = self.empty()
        with utils.open_profile(filename, 'w') as profile_handle:
            profile.save(profile_handle)

        with utils.open_profile(filename, 'r') as profile_handle:
            profile = klib.kMer.from_file(profile_handle)

        utils.test_profile(profile, counts, 4)

    def test_profile_merge(self):
        counts_left = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_right = utils.counts(utils.SEQUENCES_RIGHT, 8)

        profile_left = klib.kMer(utils.as_array(counts_left, 8))
        profile_right = klib.kMer(utils.as_array(counts_right, 8))

        profile_left.merge(profile_right)
        utils.test_profile(profile_left, counts_left + counts_right, 8)

    def test_profile_balance(self):
        counts = utils.counts(utils.SEQUENCES, 8)
        profile = klib.kMer(utils.as_array(counts, 8))
        profile.balance()

        counts.update(dict((str(Seq.reverse_complement(s)), c)
                           for s, c in counts.items()))
        utils.test_profile(profile, counts, 8)

    def test_profile_balance_palindrome(self):
        counts = utils.counts(['AATT'], 4)
        profile = klib.kMer(utils.as_array(counts, 4))
        profile.balance()

        counts.update(dict((str(Seq.reverse_complement(s)), c)
                           for s, c in counts.items()))
        utils.test_profile(profile, counts, 4)

    def _test_profile_split(self, sequences, length):
        counts = utils.counts(sequences, length)
        profile = klib.kMer(utils.as_array(counts, length))
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

    def test_profile_shrink(self):
        counts = utils.counts(utils.SEQUENCES, 8)
        profile = klib.kMer(utils.as_array(counts, 8))
        profile.shrink(1)

        counts = utils.Counter(dict((t, sum(counts[u] for u in counts
                                            if u.startswith(t)))
                                    for t in set(s[:-1] for s in counts)))
        utils.test_profile(profile, counts, 7)

    def test_profile_shrink_two(self):
        counts = utils.counts(utils.SEQUENCES, 8)
        profile = klib.kMer(utils.as_array(counts, 8))
        profile.shrink(2)

        counts = utils.Counter(dict((t, sum(counts[u] for u in counts
                                            if u.startswith(t)))
                                    for t in set(s[:-2] for s in counts)))
        utils.test_profile(profile, counts, 6)

    def test_profile_shrink_three(self):
        counts = utils.counts(utils.SEQUENCES, 8)
        profile = klib.kMer(utils.as_array(counts, 8))
        profile.shrink(3)

        counts = utils.Counter(dict((t, sum(counts[u] for u in counts
                                            if u.startswith(t)))
                                    for t in set(s[:-3] for s in counts)))
        utils.test_profile(profile, counts, 5)

    def test_profile_shrink_max(self):
        counts = utils.counts(utils.SEQUENCES, 4)
        profile = klib.kMer(utils.as_array(counts, 4))
        profile.shrink(3)

        counts = utils.Counter(dict((t, sum(counts[u] for u in counts
                                            if u.startswith(t)))
                                    for t in set(s[:-3] for s in counts)))
        utils.test_profile(profile, counts, 1)

    def test_profile_shrink_invalid(self):
        counts = utils.counts(utils.SEQUENCES, 4)
        profile = klib.kMer(utils.as_array(counts, 4))
        with pytest.raises(ValueError):
            profile.shrink(4)
