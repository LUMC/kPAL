"""
Utilities for kMer unit tests.
"""


import collections
import itertools
import os
import shutil
import tempfile

try:
    # Python 2.7 and up.
    from collections import Counter
except ImportError:
    from _counter import Counter


def fasta(sequences):
    """
    A serialization of `sequences` in FASTA format.
    """
    names = ('>sequence_%d' % i for i in itertools.count(1))
    return '\n'.join('\n'.join(entry) for entry in
                     itertools.izip(names, sequences)) + '\n'


def profile(counts, k):
    """
    A serialization of `counts` to profile format.
    """
    content = '%d\n%d\n%d\n' % (k, sum(counts.values()), len(counts))
    content += '\n'.join(str(counts[''.join(s)]) for s in
                         itertools.product('ACGT', repeat=k)) + '\n'
    return content


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
    return profile.count[count_index(sequence)]


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


class TestEnvironment():
    """
    Test class providing an isolated test environment for each test.
    """
    def setup(self):
        self.temp_dir = tempfile.mkdtemp(prefix='kmer-tests-')

    def teardown(self):
        shutil.rmtree(self.temp_dir)

    def touch(self, content=None):
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
