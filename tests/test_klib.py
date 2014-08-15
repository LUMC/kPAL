"""
Tests for the `k_mer.klib` module.
"""


from Bio import Seq

from k_mer import klib

import utils


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


class TestKlib(utils.TestEnvironment):
    def _test_profile(self, profile, counts, k):
        """
        Validate the given kMer profile, using `counts` as a reference.
        """
        assert profile.length == k
        assert profile.total == sum(counts.values())
        assert profile.non_zero == len(counts)
        for s in counts:
            assert utils.count(profile, s) == counts[s]

    def _test_analyse(self, sequences, k):
        """
        Test `k_mer.klib.analyse` with `sequences` and `k=k`.
        """
        counts = utils.counts(sequences, k)
        profile = klib.kMer()
        with open(self.touch(utils.fasta(sequences))) as fasta_handle:
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
        with open(self.touch(utils.fasta(LENGTH_60))) as fasta_handle:
            profile.analyse(fasta_handle, 4)

        filename = self.touch()
        with open(filename, 'w') as profile_handle:
            profile.save(profile_handle)

        counts = utils.counts(LENGTH_60, 4)
        with open(filename) as profile_handle:
            assert int(next(profile_handle)) == 4
            assert int(next(profile_handle)) == sum(counts.values())
            assert int(next(profile_handle)) == len(counts)
            profile_counts = [int(l) for l in profile_handle]
            for s in counts:
                assert profile_counts[utils.count_index(s)] == counts[s]

    def test_profile_load(self):
        counts = utils.counts(LENGTH_60, 4)
        filename = self.touch(utils.profile(counts, 4))

        profile = klib.kMer()
        with open(filename) as profile_handle:
            profile.load(profile_handle)
        assert profile.length == 4
        assert profile.total == sum(counts.values())
        assert profile.non_zero == len(counts)
        for s in counts:
            assert utils.count(profile, s) == counts[s]

    def test_profile_load_save(self):
        counts = utils.counts(LENGTH_60, 4)
        filename_load = self.touch(utils.profile(counts, 4))
        filename_save = self.touch()

        profile = klib.kMer()
        with open(filename_load) as profile_handle:
            profile.load(profile_handle)

        with open(filename_save, 'w') as profile_handle:
            profile.save(profile_handle)

        with open(filename_save) as f:
            assert f.read() == utils.profile(counts, 4)

    def test_profile_save_load(self):
        profile_save = klib.kMer()
        with open(self.touch(utils.fasta(LENGTH_60))) as fasta_handle:
            profile_save.analyse(fasta_handle, 4)

        filename = self.touch()
        with open(filename, 'w') as profile_handle:
            profile_save.save(profile_handle)

        profile_load = klib.kMer()
        with open(filename) as profile_handle:
            profile_load.load(profile_handle)

        assert profile_save.length == profile_load.length
        assert profile_save.total == profile_load.total
        assert profile_save.non_zero == profile_load.non_zero
        assert profile_save.count == profile_load.count

    def test_profile_balance(self):
        profile = klib.kMer()
        with open(self.touch(utils.fasta(LENGTH_60))) as fasta_handle:
            profile.analyse(fasta_handle, 8)
        profile.balance()

        counts = utils.counts(LENGTH_60, 8)
        counts.update(dict((str(Seq.reverse_complement(s)), c) for s, c in counts.items()))
        self._test_profile(profile, counts, 8)

    def test_profile_balance_palindrome(self):
        profile = klib.kMer()
        with open(self.touch(utils.fasta(['AATT']))) as fasta_handle:
            profile.analyse(fasta_handle, 4)
        profile.balance()

        counts = utils.counts(['AATT'], 4)
        counts.update(dict((str(Seq.reverse_complement(s)), c) for s, c in counts.items()))
        self._test_profile(profile, counts, 4)
