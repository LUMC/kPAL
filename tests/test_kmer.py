"""
Tests for the `k_mer.kmer` module.
"""


from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future import standard_library
from future.builtins import str, zip


import itertools
from io import open, StringIO

from Bio import Seq
import numpy as np

from k_mer import kmer

import utils

with standard_library.hooks():
    from collections import Counter


class TestKmer(utils.TestEnvironment):
    def test_convert(self):
        counts = utils.counts(utils.SEQUENCES, 8)
        filename = self.empty()
        with open(self.profile_old_format(counts, 8)) as handle:
            with utils.open_profile(filename, 'w') as profile_handle:
                kmer.convert([handle], profile_handle)
        utils.test_profile_file(filename, counts, 8)

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
                    kmer.merge(handle_left, handle_right, profile_handle)
        utils.test_profile_file(filename, counts_left + counts_right, 8)

    def test_balance(self):
        counts = utils.counts(utils.SEQUENCES, 8)
        filename = self.empty()

        with utils.open_profile(self.profile(counts, 8)) as input_handle:
            with utils.open_profile(filename, 'w') as output_handle:
                kmer.balance(input_handle, output_handle)
        counts.update(dict((utils.reverse_complement(s), c) for s, c in counts.items()))
        utils.test_profile_file(filename, counts, 8)

    def test_get_balance(self):
        counts = utils.counts(utils.SEQUENCES, 8)
        out = StringIO()

        with utils.open_profile(self.profile(counts, 8)) as input_handle:
            kmer.get_balance(input_handle, out)

        assert out.getvalue() == '1 0.669\n'

    def test_get_stats(self):
        counts = utils.counts(utils.SEQUENCES, 8)
        out = StringIO()

        with utils.open_profile(self.profile(counts, 8)) as input_handle:
            kmer.get_stats(input_handle, out)

        name, mean, std = out.getvalue().strip().split()
        assert name == '1'
        assert mean == '%.3f' % np.mean(utils.as_array(counts, 8))
        assert std == '%.3f' % np.std(utils.as_array(counts, 8))

    def test_distribution(self):
        counts = utils.counts(utils.SEQUENCES, 8)
        out = StringIO()

        with utils.open_profile(self.profile(counts, 8)) as input_handle:
            kmer.distribution(input_handle, out)

        counter = Counter(utils.as_array(counts, 8))
        assert out.getvalue() == '\n'.join('1 %i %i' % x
                                           for x in sorted(counter.items())) + '\n'

    def test_info(self):
        counts = utils.counts(utils.SEQUENCES, 8)
        out = StringIO()

        with utils.open_profile(self.profile(counts, 8, 'a')) as input_handle:
            kmer.info(input_handle, out)

        expected = 'File format version: 1.0.0\n'
        expected += 'Produced by: kMer unit tests\n\n'
        expected += 'Profile: a\n'
        expected += '- k-mer length: 8 (%d k-mers)\n' % (4**8)
        expected += '- Zero counts: %i\n' % (4**8 - len(counts))
        expected += '- Non-zero counts: %i\n' % len(counts)
        expected += '- Sum of counts: %i\n' % sum(counts.values())
        expected += '- Mean of counts: %.3f\n' % np.mean([0] * (4**8 - len(counts)) + list(counts.values()))
        expected += '- Median of counts: %i\n' % np.median([0] * (4**8 - len(counts)) + list(counts.values()))
        expected += '- Standard deviation of counts: %.3f\n' % np.std([0] * (4**8 - len(counts)) + list(counts.values()))

        assert out.getvalue() == expected

    def test_get_count(self):
        counts = utils.counts(utils.SEQUENCES, 8)
        word, count = counts.most_common(1)[0]
        out = StringIO()

        with utils.open_profile(self.profile(counts, 8, 'a')) as input_handle:
            kmer.get_count(input_handle, out, word)

        assert out.getvalue() == 'a %d\n' % count

    def test_positive(self):
        counts_left = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_right = utils.counts(utils.SEQUENCES_RIGHT, 8)
        filename_left = self.empty()
        filename_right = self.empty()

        with utils.open_profile(self.profile(counts_left, 8)) as handle_left:
            with utils.open_profile(self.profile(counts_right, 8)) as handle_right:
                with utils.open_profile(filename_left, 'w') as out_left:
                    with utils.open_profile(filename_right, 'w') as out_right:
                        kmer.positive(handle_left, handle_right, out_left, out_right)

        utils.test_profile_file(filename_left, Counter(s for s in counts_left.elements()
                                                             if s in counts_right), 8)
        utils.test_profile_file(filename_right, Counter(s for s in counts_right.elements()
                                                              if s in counts_left), 8)

    def test_scale(self):
        counts_left = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_right = utils.counts(utils.SEQUENCES_RIGHT, 8)
        filename_left = self.empty()
        filename_right = self.empty()

        with utils.open_profile(self.profile(counts_left, 8)) as handle_left:
            with utils.open_profile(self.profile(counts_right, 8)) as handle_right:
                with utils.open_profile(filename_left, 'w') as out_left:
                    with utils.open_profile(filename_right, 'w') as out_right:
                        kmer.scale(handle_left, handle_right, out_left, out_right)

        if sum(counts_left.values()) < sum(counts_right.values()):
            scale_left = sum(counts_right.values()) / sum(counts_left.values())
            scale_right = 1.0
        else:
            scale_left = 1.0
            scale_right = sum(counts_left.values()) / sum(counts_right.values())

        for s in counts_left:
            counts_left[s] *= scale_left
        for s in counts_right:
            counts_right[s] *= scale_right

        utils.test_profile_file(filename_left, counts_left, 8)
        utils.test_profile_file(filename_right, counts_right, 8)

    def test_shrink(self):
        counts = utils.counts(utils.SEQUENCES, 8)
        filename = self.empty()

        with utils.open_profile(self.profile(counts, 8)) as input_handle:
            with utils.open_profile(filename, 'w') as output_handle:
                kmer.shrink(input_handle, output_handle, 1)

        counts = Counter(dict((t, sum(counts[u] for u in counts
                                            if u.startswith(t)))
                                    for t in set(s[:-1] for s in counts)))
        utils.test_profile_file(filename, counts, 7)

    def test_shuffle(self):
        # See test_klib.profile_shuffle
        counts = utils.counts(utils.SEQUENCES, 2)
        filename = self.empty()

        with utils.open_profile(self.profile(counts, 2)) as input_handle:
            with utils.open_profile(filename, 'w') as output_handle:
                np.random.seed(100)
                kmer.shuffle(input_handle, output_handle)

        counts = dict(zip([''.join(s) for s in itertools.product('ACGT', repeat=2)],
                          [13,  7,  6, 18, 12,  1, 13, 17, 16, 12, 23, 27, 24, 17, 18, 12]))
        utils.test_profile_file(filename, counts, 2)

    def test_smooth(self):
        # See test_kdifflib.test_kmerdiff_dynamic_smooth
        counts_left = Counter(['AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TG', 'TT'])
        counts_right = Counter(['AC', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT'])
        filename_left = self.empty()
        filename_right = self.empty()

        with utils.open_profile(self.profile(counts_left, 2)) as handle_left:
            with utils.open_profile(self.profile(counts_right, 2)) as handle_right:
                with utils.open_profile(filename_left, 'w') as out_left:
                    with utils.open_profile(filename_right, 'w') as out_right:
                        kmer.smooth(handle_left, handle_right, out_left, out_right, summary='min')

        counts_left = Counter(['AA', 'AA', 'AA', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TA', 'TA'])
        counts_right = Counter(['AA', 'AA', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TA', 'TA', 'TA'])

        utils.test_profile_file(filename_left, counts_left, 2)
        utils.test_profile_file(filename_right, counts_right, 2)

    def test_pair_diff(self):
        counts_left = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_right = utils.counts(utils.SEQUENCES_RIGHT, 8)
        out = StringIO()

        with utils.open_profile(self.profile(counts_left, 8, 'left')) as handle_left:
            with utils.open_profile(self.profile(counts_right, 8, 'right')) as handle_right:
                kmer.pair_diff(handle_left, handle_right, out)

        assert out.getvalue() == 'left right %.3f\n' % 0.4626209322

    def test_matrix_diff(self):
        counts_left = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_right = utils.counts(utils.SEQUENCES_RIGHT, 8)
        out = StringIO()

        with utils.open_profile(self.multi_profile(8,
                                                   [counts_left,
                                                    counts_right,
                                                    counts_left],
                                                   ['a', 'b', 'c'])) as handle:
                    kmer.matrix_diff(handle, out)

        assert out.getvalue().strip().split('\n') == ['3', 'a', 'b', 'c', '0.463', '0.000 0.463']
