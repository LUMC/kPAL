"""
k-mer base library.
"""

import argparse
import itertools
import math
import os
import random
import re
import sys

from Bio import Seq, SeqIO
import numpy as np

from . import metrics

class Profile():
    """
    Handle k-mer counts.
    """
    __nucleotide_to_binary = {
        'A': 0x00, 'a': 0x00,
        'C': 0x01, 'c': 0x01,
        'G': 0x02, 'g': 0x02,
        'T': 0x03, 't': 0x03
    }
    """ Conversion table form nucleotide to binary. """
    __binary_to_nucleotide = {
        0x00: 'A',
        0x01: 'C',
        0x02: 'G',
        0x03: 'T'
    }
    """ Conversion table form binary to nucleotide. """

    @classmethod
    def from_file(cls, handle, name=None):
        """
        Load the k-mer profile from a file.

        :arg handle: Open readable k-mer profile file handle.
        :type handle: h5py.File
        :arg name: Profile name.
        :arg name: str

        :return: A k-mer profile.
        :rtype: Profile
        """
        name = name or sorted(handle['profiles'].keys())[0]
        counts = handle['profiles/%s' % name][:]
        return cls(counts, name=name)
    #from_file

    @classmethod
    def from_file_old_format(cls, handle, name=None):
        """
        Load the k-mer profile from a file in the old plaintext format.

        :arg handle: Open readable k-mer profile file handle (old format).
        :type handle: file-like object
        :arg name: Profile name.
        :arg name: str

        :return: A k-mer profile.
        :rtype: Profile
        """
        # Ignore lines with length, total, nonzero.
        for _ in range(3):
            next(handle)

        counts = np.loadtxt(handle, dtype='int64')
        return cls(counts, name=name)
    #from_file_old_format

    @classmethod
    def from_fasta(cls, handle, length, name=None):
        """
        Create a k-mer profile from a FASTA file by counting all k-mers in
        each line.

        :arg handle: Open readable FASTA file handle.
        :type handle: file-like object
        :arg length: Length of the k-mers.
        :type length: int
        :arg name: Profile name.
        :arg name: str

        :return: A k-mer profile.
        :rtype: Profile
        """
        number = 4 ** length
        bitmask = number - 0x01
        counts = [0] * number
        alphabet = re.compile("[^%s]" %
            ''.join(cls.__nucleotide_to_binary.keys()))

        for record in SeqIO.parse(handle, "fasta"):
            for sequence in alphabet.split(str(record.seq)):
                if len(sequence) >= length:
                    binary = 0x00

                    # Calculate the binary representation of a k-mer.
                    for i in sequence[:length]:
                        binary = ((binary << 2) |
                            cls.__nucleotide_to_binary[i])
                    counts[binary] += 1

                    # Calculate the binary representation of the next k-mer.
                    for i in sequence[length:]:
                        binary = ((binary << 2) |
                            cls.__nucleotide_to_binary[i]) & bitmask
                        counts[binary] += 1

        return cls(np.array(counts, dtype='int64'), name=name)
    #from_fasta

    def __init__(self, counts, name=None):
        """
        Low-level constructor for k-mer profile objects. Users should
        generally use one of the profile construction methods:

        - `kMer.from_file`
        - `kMer.from_file_old_format`
        - `kMer.from_fasta`

        :arg counts: NumPy array of integers where each element is the count
          for a k-mer. Ordering is alphabetically by the k-mer.
        :type counts: np.ndarray(int)
        :arg name: Profile name.
        :type name: str
        """
        self.length = int(math.log(len(counts), 4))
        self.counts = counts
        self.name = name
    #__init__

    @property
    def number(self):
        """
        Number of possible k-mers with this length.
        """
        return len(self.counts)

    @property
    def non_zero(self):
        """
        Number k-mers with a non-zero count.
        """
        return np.count_nonzero(self.counts)

    @property
    def total(self):
        """
        Sum of k-mer counts.
        """
        return self.counts.sum()

    @property
    def mean(self):
        """
        Mean of k-mer counts.
        """
        return self.counts.mean()

    @property
    def median(self):
        """
        Median of k-mer counts.
        """
        return np.median(self.counts)

    @property
    def std(self):
        """
        Standard deviation of k-mer counts.
        """
        return self.counts.std()

    def save(self, handle, name=None):
        """
        Save the k-mer counts to a file.

        :arg handle: Open writeable k-mer profile file handle.
        :type handle: h5py.File
        :arg name: Profile name in the file. If not provided, the current
          profile name is used, or the first available number from 1
          consecutively if the profile has no name.
        :type name: str

        :return: Profile name in the file.
        :rtype: str
        """
        name = name or self.name or next(str(n) for n in itertools.count(1)
                                         if str(n) not in handle['profiles'])

        profile = handle.create_dataset('profiles/%s' % name, data=self.counts,
                                        dtype='int64', compression='gzip')
        profile.attrs['length'] = self.length
        profile.attrs['total'] = self.total
        profile.attrs['non_zero'] = self.non_zero
        profile.attrs['mean'] = self.mean
        profile.attrs['median'] = self.median
        profile.attrs['std'] = self.std
        handle.flush()

        return name
    #save

    def merge(self, profile, merger=metrics.mergers["sum"]):
        """
        Merge two profiles.

        :arg profile: Another k-mer profile.
        :type profile: Profile
        :arg merger: Merge function.
        :type merger: function
        """
        self.counts = merger(self.counts, profile.counts)
    #merge

    def balance(self):
        """
        Add the counts of the reverse complement of a k-mer to the k-mer and
        vice versa.
        """
        for i in range(self.number):
            i_rc = self.reverse_complement(i)

            if i < i_rc:
                temp = self.counts[i]
                self.counts[i] += self.counts[i_rc]
                self.counts[i_rc] += temp
            elif i == i_rc:
                self.counts[i] += self.counts[i]
            #if
        #for
    #balance

    def split(self):
        """
        Split the profile into two lists, every position in the first list has
        its reverse complement in the same position in the second list and vice
        versa. All counts are doubled, so we can equaly distribute palindrome
        counts over both lists.

        Note that the returned counts are not k-mer profiles. They can be used
        to show the balance of the original profile by calculating the distance
        between them.

        :return: The doubled forward and reverse complement counts.
        :rtype: np.ndarray(int), np.ndarray(int)
        """
        forward = []
        reverse = []

        for i in range(self.number):
            i_rc = self.reverse_complement(i)

            if i < i_rc:
                forward.append(self.counts[i] * 2)
                reverse.append(self.counts[i_rc] * 2)
            elif i == i_rc:
                forward.append(self.counts[i])
                reverse.append(self.counts[i])
            #if
        #for

        return np.array(forward), np.array(reverse)
    #split

    def shrink(self, factor=1):
        """
        Shrink the profile, effectively reducing the value of k.

        Note that this operation may give slightly different values than
        indexing on a lower k directly.

        :arg factor: Shrinking factor.
        :type factor: int
        """
        if self.length <= factor:
            raise ValueError(
                "Reduction factor should be smaller than k-mer size.")

        # I don't know how to do this operation directly with NumPy
        # instructions without creating a complete copy of the counts in
        # memory (e.g., using `np.split(self.counts, ...)`).
        # So instead we just build the new counts one at a time. Note that
        # the builtin Python `sum` is actually faster than `np.sum` here,
        # since the arrays are very small.
        merge_size = 4 ** factor
        new_counts = (sum(self.counts[i:i + merge_size])
                      for i in range(0, self.number, merge_size))
        self.counts = np.fromiter(new_counts, dtype='int64')
        self.length -= factor
    #shrink

    def shuffle(self):
        """
        Randomise the profile.
        """
        np.random.shuffle(self.counts)
    #shuffle

    def dna_to_binary(self, sequence):
        """
        Convert a string of DNA to an integer.

        :arg sequence: DNA sequence.
        :type sequence: str

        :return: Binary representation of `sequence`.
        :rtype: int
        """
        result = 0x00

        for i in sequence:
            result <<= 2
            result |= self.__nucleotide_to_binary[i]
        #for

        return result
    #dna_to_binary

    def binary_to_dna(self, number):
        """
        Convert an integer to a DNA string.

        :arg number: Binary representation of a DNA sequence.
        :type number: int

        :returns DNA string corresponding to `number`.
        :rtype: str
        """
        sequence = ""

        for i in range(self.length):
            sequence += self.__binary_to_nucleotide[number & 0x03]
            number >>= 2
        #while

        return sequence[::-1]
    #binary_to_dna

    def reverse_complement(self, number):
        """
        Calculate the reverse complement of a DNA sequence in a binary
        representation.

        :arg number: Binary representation of a DNA sequence.
        :type number: int

        :return: Binary representation of the reverse complement of the
          sequence corresponding to `number`.
        :rtype: int
        """
        number = ~number
        result = 0x00

        for i in range(self.length):
            result = (result << 2) | (number & 0x03)
            number >>= 2
        #for

        return result
    #reverse_complement

    def _ratios_matrix(self):
        """
        Calculate all relative frequencies of k-mers. If a division by 0
        occurs, the frequency will be set to -1.0.

        :return: A matrix with relative frequencies.
        :rtype: float[][]
        """
        # Perhaps we can use something like this for a future normalization
        # operation.
        # Todo: Do this directly with NumPy.
        ratios = []
        for i in range(self.number):
            ratios.append(self.number * [0.0])

        # Fill the matrix.
        for i in range(self.number):
            for j in range(self.number):
                if self.counts[j]:
                    ratios[i][j] = (float(self.counts[i]) /
                        self.counts[j]) / self.total
                else:
                    ratios[i][j] = -1.0
            #for

        return ratios
    #ratios_matrix

    def _freq_diff_matrix(self):
        """
        Calculate all frequency differences of k-mers.

        :return: A matrix with frequency differences.
        :rtype: float[][]
        """
        ratios = []
        for i in range(self.number):
            ratios.append(self.number * [0])

        # Fill the matrix.
        for i in range(self.number):
            for j in range(self.number):
                if self.counts[j]:
                    ratios[i][j] = float(abs(self.counts[i] -
                        self.counts[j])) / self.total

        return ratios
    #freq_diff_matrix

    def print_counts(self):
        """
        Print the k-mer counts.
        """
        for i in range(self.number):
            print self.binary_to_dna(i), self.counts[i]
    #print_counts

    def _print_ratios(self, ratios):
        """
        Print a ratios matrix.

        :arg ratios: A matrix with relative frequencies.
        :type ratios: float[][]
        """
        # The header.
        print (self.length + 1) * ' ',
        for i in range(self.number):
            print "%s%s" % (self.binary_to_dna(i), 2 * ' '),
        print

        # The matrix.
        for i in range(self.number):
            print "%s" % self.binary_to_dna(i),
            for j in range(self.number):
                print ("%%.%if" % self.length) % ratios[i][j],
            print
        #for
    #print_ratios
#Profile
