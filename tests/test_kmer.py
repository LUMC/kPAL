"""
Tests for the `kpal.kmer` module.
"""


from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future import standard_library
from future.builtins import str, zip


import itertools
from io import open, StringIO

from Bio import Seq
import numpy as np

from kpal import kmer

import utils

with standard_library.hooks():
    from collections import Counter


class TestKmer(utils.TestEnvironment):
    def test_main_info(self, capsys):
        # For the `capsys` fixture, see:
        # http://pytest.org/latest/capture.html

        counts = utils.counts(utils.SEQUENCES, 8)
        filename = self.profile(counts, 8, 'a')

        kmer.main(['info', filename])

        out, err = capsys.readouterr()

        expected = 'File format version: 1.0.0\n'
        expected += 'Produced by: kMer unit tests\n\n'
        expected += 'Profile: a\n'
        expected += '- k-mer length: 8 (%d k-mers)\n' % (4**8)
        expected += '- Zero counts: %i\n' % (4**8 - len(counts))
        expected += '- Non-zero counts: %i\n' % len(counts)
        expected += '- Sum of counts: %i\n' % sum(counts.values())
        expected += '- Mean of counts: %.3f\n' % np.mean([0] * (4**8 - len(counts)) + list(counts.values()))
        expected += '- Median of counts: %.3f\n' % np.median([0] * (4**8 - len(counts)) + list(counts.values()))
        expected += '- Standard deviation of counts: %.3f\n' % np.std([0] * (4**8 - len(counts)) + list(counts.values()))

        assert out == expected

    def test_convert(self):
        counts = utils.counts(utils.SEQUENCES, 8)
        filename = self.empty()
        with open(self.profile_old_format(counts, 8)) as handle:
            with utils.open_profile(filename, 'w') as profile_handle:
                kmer.convert([handle], profile_handle)
        utils.test_profile_file(filename, counts, 8)

    def test_cat(self):
        counts_a = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_b = utils.counts(utils.SEQUENCES_RIGHT, 8)
        filename = self.empty()

        with utils.open_profile(self.profile(counts_a, 8, name='a')) as handle_a:
            with utils.open_profile(self.profile(counts_b, 8, name='b')) as handle_b:
                with utils.open_profile(filename, 'w') as profile_handle:
                    kmer.cat([handle_a, handle_b], profile_handle)
        utils.test_profile_file(filename, counts_a, 8, name='a')
        utils.test_profile_file(filename, counts_b, 8, name='b')

    def test_cat_prefixes(self):
        counts_a = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_b = utils.counts(utils.SEQUENCES_RIGHT, 8)
        filename = self.empty()

        with utils.open_profile(self.profile(counts_a, 8, name='X')) as handle_a:
            with utils.open_profile(self.profile(counts_b, 8, name='X')) as handle_b:
                with utils.open_profile(filename, 'w') as profile_handle:
                    kmer.cat([handle_a, handle_b], profile_handle, prefixes=['a_', 'b_'])
        utils.test_profile_file(filename, counts_a, 8, name='a_X')
        utils.test_profile_file(filename, counts_b, 8, name='b_X')

    def test_count(self):
        counts = utils.counts(utils.SEQUENCES, 8)
        filename = self.empty()
        with open(self.fasta(utils.SEQUENCES)) as fasta_handle:
            with utils.open_profile(filename, 'w') as profile_handle:
                kmer.count([fasta_handle], profile_handle, 8)
        utils.test_profile_file(filename, counts, 8)

    def test_count_multi(self):
        counts_left = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_right = utils.counts(utils.SEQUENCES_RIGHT, 8)
        filename = self.empty()
        with open(self.fasta(utils.SEQUENCES_LEFT)) as handle_left:
            with open(self.fasta(utils.SEQUENCES_RIGHT)) as handle_right:
                with utils.open_profile(filename, 'w') as profile_handle:
                    kmer.count([handle_left, handle_right], profile_handle, 8, names=['a', 'b'])
        utils.test_profile_file(filename, counts_left, 8, name='a')
        utils.test_profile_file(filename, counts_right, 8, name='b')

    def test_merge(self):
        counts_left = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_right = utils.counts(utils.SEQUENCES_RIGHT, 8)
        filename = self.empty()

        with utils.open_profile(self.profile(counts_left, 8)) as handle_left:
            with utils.open_profile(self.profile(counts_right, 8)) as handle_right:
                with utils.open_profile(filename, 'w') as profile_handle:
                    kmer.merge(handle_left, handle_right, profile_handle)
        utils.test_profile_file(filename, counts_left + counts_right, 8)

    def test_merge_xor(self):
        counts_left = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_right = utils.counts(utils.SEQUENCES_RIGHT, 8)
        filename = self.empty()

        with utils.open_profile(self.profile(counts_left, 8)) as handle_left:
            with utils.open_profile(self.profile(counts_right, 8)) as handle_right:
                with utils.open_profile(filename, 'w') as profile_handle:
                    kmer.merge(handle_left, handle_right, profile_handle, merger='xor')

        counts_xor = counts_left + counts_right
        for s in set(counts_left) & set(counts_right):
            del counts_xor[s]

        utils.test_profile_file(filename, counts_xor, 8)

    def test_merge_custom_expr(self):
        counts_left = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_right = utils.counts(utils.SEQUENCES_RIGHT, 8)
        filename = self.empty()

        with utils.open_profile(self.profile(counts_left, 8)) as handle_left:
            with utils.open_profile(self.profile(counts_right, 8)) as handle_right:
                with utils.open_profile(filename, 'w') as profile_handle:
                    kmer.merge(handle_left, handle_right, profile_handle, custom_merger='(left + right) * np.logical_xor(left, right)')

        counts_xor = counts_left + counts_right
        for s in set(counts_left) & set(counts_right):
            del counts_xor[s]

        utils.test_profile_file(filename, counts_xor, 8)

    def test_merge_custom_name(self):
        counts_left = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_right = utils.counts(utils.SEQUENCES_RIGHT, 8)
        filename = self.empty()

        with utils.open_profile(self.profile(counts_left, 8)) as handle_left:
            with utils.open_profile(self.profile(counts_right, 8)) as handle_right:
                with utils.open_profile(filename, 'w') as profile_handle:
                    kmer.merge(handle_left, handle_right, profile_handle, custom_merger='numpy.multiply')

        counts_mult = Counter(dict((s, counts_left[s] * counts_right[s])
                                   for s in set(counts_left) & set(counts_right)))

        utils.test_profile_file(filename, counts_mult, 8)

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
            kmer.get_balance(input_handle, out, precision=3)

        assert out.getvalue() == '1 0.669\n'

    def test_get_stats(self):
        counts = utils.counts(utils.SEQUENCES, 8)
        out = StringIO()

        with utils.open_profile(self.profile(counts, 8)) as input_handle:
            kmer.get_stats(input_handle, out)

        name, mean, std = out.getvalue().strip().split()
        assert name == '1'
        assert mean == '%.10f' % np.mean(utils.as_array(counts, 8))
        assert std == '%.10f' % np.std(utils.as_array(counts, 8))

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
        expected += '- Median of counts: %.3f\n' % np.median([0] * (4**8 - len(counts)) + list(counts.values()))
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
        # See test_kdistlib.test_ProfileDistance_dynamic_smooth
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

    def test_smooth_custom_expr(self):
        # See test_kdistlib.test_ProfileDistance_dynamic_smooth
        counts_left = Counter(['AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TG', 'TT'])
        counts_right = Counter(['AC', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT'])
        filename_left = self.empty()
        filename_right = self.empty()

        with utils.open_profile(self.profile(counts_left, 2)) as handle_left:
            with utils.open_profile(self.profile(counts_right, 2)) as handle_right:
                with utils.open_profile(filename_left, 'w') as out_left:
                    with utils.open_profile(filename_right, 'w') as out_right:
                        kmer.smooth(handle_left, handle_right, out_left, out_right, custom_summary='np.max(values)')

    def test_smooth_custom_name(self):
        # See test_kdistlib.test_ProfileDistance_dynamic_smooth
        counts_left = Counter(['AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TG', 'TT'])
        counts_right = Counter(['AC', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT'])
        filename_left = self.empty()
        filename_right = self.empty()

        with utils.open_profile(self.profile(counts_left, 2)) as handle_left:
            with utils.open_profile(self.profile(counts_right, 2)) as handle_right:
                with utils.open_profile(filename_left, 'w') as out_left:
                    with utils.open_profile(filename_right, 'w') as out_right:
                        kmer.smooth(handle_left, handle_right, out_left, out_right, custom_summary='numpy.max')

    def test_distance(self):
        counts_left = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_right = utils.counts(utils.SEQUENCES_RIGHT, 8)
        out = StringIO()

        with utils.open_profile(self.profile(counts_left, 8, 'left')) as handle_left:
            with utils.open_profile(self.profile(counts_right, 8, 'right')) as handle_right:
                kmer.distance(handle_left, handle_right, out)

        assert out.getvalue() == 'left right %.10f\n' % 0.4626209323

    def test_distance_smooth(self):
        counts_left = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_right = utils.counts(utils.SEQUENCES_RIGHT, 8)
        out = StringIO()

        with utils.open_profile(self.profile(counts_left, 8, 'left')) as handle_left:
            with utils.open_profile(self.profile(counts_right, 8, 'right')) as handle_right:
                kmer.distance(handle_left, handle_right, out, do_smooth=True, precision=3)

        assert out.getvalue() == 'left right 0.077\n'

    def test_distance_smooth_average(self):
        counts_left = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_right = utils.counts(utils.SEQUENCES_RIGHT, 8)
        out = StringIO()

        with utils.open_profile(self.profile(counts_left, 8, 'left')) as handle_left:
            with utils.open_profile(self.profile(counts_right, 8, 'right')) as handle_right:
                kmer.distance(handle_left, handle_right, out, do_smooth=True,
                              precision=3, summary='average')

        assert out.getvalue() == 'left right 0.474\n'

    def test_distance_smooth_expr(self):
        counts_left = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_right = utils.counts(utils.SEQUENCES_RIGHT, 8)
        out = StringIO()

        with utils.open_profile(self.profile(counts_left, 8, 'left')) as handle_left:
            with utils.open_profile(self.profile(counts_right, 8, 'right')) as handle_right:
                kmer.distance(handle_left, handle_right, out, do_smooth=True,
                              precision=3, custom_summary='np.max(values)')

        assert out.getvalue() == 'left right 0.474\n'

    def test_distance_smooth_name(self):
        counts_left = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_right = utils.counts(utils.SEQUENCES_RIGHT, 8)
        out = StringIO()

        with utils.open_profile(self.profile(counts_left, 8, 'left')) as handle_left:
            with utils.open_profile(self.profile(counts_right, 8, 'right')) as handle_right:
                kmer.distance(handle_left, handle_right, out, do_smooth=True,
                              precision=3, custom_summary='numpy.max')

        assert out.getvalue() == 'left right 0.474\n'

    def test_distance_pairwise_expr(self):
        counts_left = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_right = utils.counts(utils.SEQUENCES_RIGHT, 8)
        out = StringIO()

        with utils.open_profile(self.profile(counts_left, 8, 'left')) as handle_left:
            with utils.open_profile(self.profile(counts_right, 8, 'right')) as handle_right:
                kmer.distance(handle_left, handle_right, out, precision=3,
                              custom_pairwise='abs(left - right) / (left + right + 1000)')

        assert out.getvalue() == 'left right 0.001\n'

    def test_distance_pairwise_name(self):
        counts_left = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_right = utils.counts(utils.SEQUENCES_RIGHT, 8)
        out = StringIO()

        with utils.open_profile(self.profile(counts_left, 8, 'left')) as handle_left:
            with utils.open_profile(self.profile(counts_right, 8, 'right')) as handle_right:
                kmer.distance(handle_left, handle_right, out, precision=3,
                              custom_pairwise='numpy.multiply')

        assert out.getvalue() == 'left right 0.084\n'

    def test_distance_matrix(self):
        counts_left = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_right = utils.counts(utils.SEQUENCES_RIGHT, 8)
        out = StringIO()

        with utils.open_profile(self.multi_profile(8,
                                                   [counts_left,
                                                    counts_right,
                                                    counts_left],
                                                   ['a', 'b', 'c'])) as handle:
                    kmer.distance_matrix(handle, out, precision=3)

        assert out.getvalue().strip().split('\n') == ['3', 'a', 'b', 'c', '0.463', '0.000 0.463']

    def test_distance_matrix_smooth(self):
        counts_left = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_right = utils.counts(utils.SEQUENCES_RIGHT, 8)
        out = StringIO()

        with utils.open_profile(self.multi_profile(8,
                                                   [counts_left,
                                                    counts_right,
                                                    counts_left],
                                                   ['a', 'b', 'c'])) as handle:
                    kmer.distance_matrix(handle, out, do_smooth=True, precision=3)

        assert out.getvalue().strip().split('\n') == ['3', 'a', 'b', 'c', '0.077', '0.000 0.077']

    def test_distance_matrix_smooth_average(self):
        counts_left = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_right = utils.counts(utils.SEQUENCES_RIGHT, 8)
        out = StringIO()

        with utils.open_profile(self.multi_profile(8,
                                                   [counts_left,
                                                    counts_right,
                                                    counts_left],
                                                   ['a', 'b', 'c'])) as handle:
                    kmer.distance_matrix(handle, out, do_smooth=True,
                                         summary='average', precision=3)

        assert out.getvalue().strip().split('\n') == ['3', 'a', 'b', 'c', '0.474', '0.000 0.474']

    def test_distance_matrix_smooth_expr(self):
        counts_left = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_right = utils.counts(utils.SEQUENCES_RIGHT, 8)
        out = StringIO()

        with utils.open_profile(self.multi_profile(8,
                                                   [counts_left,
                                                    counts_right,
                                                    counts_left],
                                                   ['a', 'b', 'c'])) as handle:
                    kmer.distance_matrix(handle, out, do_smooth=True, precision=3,
                                         custom_summary='np.max(values)')

        assert out.getvalue().strip().split('\n') == ['3', 'a', 'b', 'c', '0.474', '0.000 0.474']

    def test_distance_matrix_smooth_name(self):
        counts_left = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_right = utils.counts(utils.SEQUENCES_RIGHT, 8)
        out = StringIO()

        with utils.open_profile(self.multi_profile(8,
                                                   [counts_left,
                                                    counts_right,
                                                    counts_left],
                                                   ['a', 'b', 'c'])) as handle:
                    kmer.distance_matrix(handle, out, do_smooth=True, precision=3,
                                         custom_summary='numpy.max')

        assert out.getvalue().strip().split('\n') == ['3', 'a', 'b', 'c', '0.474', '0.000 0.474']

    def test_distance_matrix_pairwise_expr(self):
        counts_left = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_right = utils.counts(utils.SEQUENCES_RIGHT, 8)
        out = StringIO()

        with utils.open_profile(self.multi_profile(8,
                                                   [counts_left,
                                                    counts_right,
                                                    counts_left],
                                                   ['a', 'b', 'c'])) as handle:
                    kmer.distance_matrix(handle, out, precision=3,
                                         custom_pairwise='abs(left - right) / (left + right + 1000)')

        assert out.getvalue().strip().split('\n') == ['3', 'a', 'b', 'c', '0.001', '0.000 0.001']

    def test_distance_matrix_pairwise_name(self):
        counts_left = utils.counts(utils.SEQUENCES_LEFT, 8)
        counts_right = utils.counts(utils.SEQUENCES_RIGHT, 8)
        out = StringIO()

        with utils.open_profile(self.multi_profile(8,
                                                   [counts_left,
                                                    counts_right,
                                                    counts_left],
                                                   ['a', 'b', 'c'])) as handle:
                    kmer.distance_matrix(handle, out, precision=3,
                                         custom_pairwise='numpy.multiply')

        assert out.getvalue().strip().split('\n') == ['3', 'a', 'b', 'c', '0.084', '1.206 0.084']
