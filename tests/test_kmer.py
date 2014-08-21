"""
Tests for the `k_mer.kmer` module.
"""


from Bio import Seq

from k_mer import kmer

import utils


class TestKmer(utils.TestEnvironment):
    def test_index(self):
        counts = utils.counts(utils.SEQUENCES, 8)
        filename = self.empty()
        with open(self.fasta(utils.SEQUENCES)) as fasta_handle:
            with utils.open_profile(filename, 'w') as profile_handle:
                kmer.index([fasta_handle], profile_handle, 8)
        utils.test_profile_file(filename, counts, 8)

    def test_merge(self):
        counts_left = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_right = utils.counts(utils.SEQUENCES_RIGHT, 8)
        filename = self.empty()

        with utils.open_profile(self.profile(counts_left, 8)) as handle_left:
            with utils.open_profile(self.profile(counts_right, 8)) as handle_right:
                with utils.open_profile(filename, 'w') as profile_handle:
                    kmer.merge([handle_left, handle_right], profile_handle)
        utils.test_profile_file(filename, counts_left + counts_right, 8)

    def test_balance(self):
        counts = utils.counts(utils.SEQUENCES, 8)
        filename = self.empty()

        with utils.open_profile(self.profile(counts, 8)) as input_handle:
            with utils.open_profile(filename, 'w') as output_handle:
                kmer.balance(input_handle, output_handle)
        counts.update(dict((str(Seq.reverse_complement(s)), c) for s, c in counts.items()))
        utils.test_profile_file(filename, counts, 8)

    def test_get_balance(self, capsys):
        # For the `capsys` fixture, see:
        # http://pytest.org/latest/capture.html
        counts = utils.counts(utils.SEQUENCES, 8)
        filename = self.empty()

        with utils.open_profile(self.profile(counts, 8)) as input_handle:
            kmer.get_balance(input_handle)

        out, err = capsys.readouterr()
        assert out == '0.669\n'

    def test_shrink(self):
        counts = utils.counts(utils.SEQUENCES, 8)
        filename = self.empty()

        with utils.open_profile(self.profile(counts, 8)) as input_handle:
            with utils.open_profile(filename, 'w') as output_handle:
                kmer.shrink(input_handle, output_handle, 1)
        counts = utils.Counter(dict((t, sum(counts[u] for u in counts
                                            if u.startswith(t)))
                                    for t in set(s[:-1] for s in counts)))
        utils.test_profile_file(filename, counts, 7)
