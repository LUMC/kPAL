"""
Tests for the `k_mer.kmer` module.
"""

from Bio import Seq

from k_mer import kmer

import utils


# Some 60-mers.
SEQUENCES_A = ['TTACAATGATTAGGTCCACAGCTCTGAGCAACGCGCAGACGTCACATACTTCAAAACCCA',
               'TAAAAACTATATAAGAAATCGAATTTTCTCTTAATGGTAGCAGCTACCGTAAAATCTATG',
               'TACGCCTATATATCTTTGACTAAGCATTTATGTATTACATACTAACCAAAATTACTGTCT',
               'TACTAAGTTTTGAAGGACAGCACATCACCTGCGCAATATCGGTGTCACCCCATAGCTCCT']

# Some 60-mers.
SEQUENCES_B = ['GTACATGATAGGTCCACAGCTCTGAGCAAGGCAGACGTCCATACTTAAAACCCAGACTGC',
               'TAAACTAAAAGAAAGAATTTTTTTAATGGTAGACTACCTAAAATTATGTCTCTTAGTCCT',
               'TATCTTTACCTATATATTTGACTAAGATTTAGTATTACTACTACCTAAAATTATGTCTCT',
               'TACTATGTCTTGAAGGACAGCACCTGACCTCCCCCTGCAAGGTGTCATCCCCAAGCTGGT']


class TestKmer(utils.TestEnvironment):
    def _test_profile(self, filename, counts, k):
        """
        Validate the given kMer profile, using `counts` as a reference.
        """
        with open(filename) as f:
            assert f.read() == utils.profile(counts, k)

    def test_index(self):
        counts = utils.counts(SEQUENCES_A, 8)
        filename = self.touch()
        with open(self.touch(utils.fasta(SEQUENCES_A))) as fasta_handle:
            with open(filename, 'w') as profile_handle:
                kmer.index(fasta_handle, profile_handle, 8)
        self._test_profile(filename, counts, 8)

    def test_merge(self):
        counts_left = utils.counts(SEQUENCES_A, 8)
        counts_right = utils.counts(SEQUENCES_B, 8)
        filename = self.touch()
        with open(self.touch(utils.profile(counts_left, 8))) as handle_left:
            with open(self.touch(utils.profile(counts_right, 8))) as handle_right:
                with open(filename, 'w') as profile_handle:
                    kmer.merge([handle_left, handle_right], profile_handle)
        self._test_profile(filename, counts_left + counts_right, 8)

    def test_balance(self):
        counts = utils.counts(SEQUENCES_A, 8)
        filename = self.touch()
        with open(self.touch(utils.profile(counts, 8))) as input_handle:
            with open(filename, 'w') as output_handle:
                kmer.balance(input_handle, output_handle)
        counts.update(dict((str(Seq.reverse_complement(s)), c) for s, c in counts.items()))
        self._test_profile(filename, counts, 8)
