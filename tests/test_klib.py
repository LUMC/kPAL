"""
Tests for the `k_mer.klib` module.
"""


from Bio import Seq

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
