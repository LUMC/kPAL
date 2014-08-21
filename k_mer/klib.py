#!/usr/bin/python

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

class kMer():
    """
    Handle the counting of k-mers in a fasta file.
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
        Load the k-mer table from a file.

        :arg handle: Open handle to a file.
        :type handle: stream
        """
        name = name or handle['profiles'].keys()[0]
        counts = handle['profiles/%s' % name][:]
        return cls(counts, name=name)
    #from_file

    @classmethod
    def from_fasta(cls, handle, length, name=None):
        """
        Read a fasta file and count all k-mers in each line.

        :arg handle: An open file handle to a fasta file.
        :type handle: stream
        :arg length: Length of the k-mers.
        :type length: int
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

        return cls(np.array(counts), name=name)
    #from_fasta

    def __init__(self, counts, name=None):
        """
        Initialise the class instance.
        """
        self.length = int(math.log(len(counts), 4))
        self.counts = counts
        self.name = name
    #__init__

    @property
    def number(self):
        return len(self.counts)

    @property
    def total(self):
        return self.counts.sum()

    @property
    def non_zero(self):
        return np.count_nonzero(self.counts)

    def save(self, handle):
        """
        Save the k-mer table in a file.

        :arg handle: Writable open handle to a file.
        :type handle: stream
        """
        name = self.name or next(str(n) for n in itertools.count(1)
                                 if str(n) not in handle['profiles'])

        profile = handle.create_dataset('profiles/%s' % name, data=self.counts,
                                        dtype='int64', compression='gzip')
        profile.attrs['length'] = self.length
        profile.attrs['total'] = self.total
        profile.attrs['non_zero'] = self.non_zero
        handle.flush()
    #save

    def merge(self, profile, merger=metrics.mergers["sum"]):
        """
        Merge two profiles.

        :arg profile: An other k-mer profile.
        :type profile: object(kMer)
        :arg merger: Merge function.
        :type merger: function
        """
        for i in range(self.number):
            self.counts[i] = merger(self.counts[i], profile.counts[i])
        #for
    #merge

    def balance(self):
        """
        Add the counts of the reverse complement of a k-mer to the k-mer and
        vice versa.
        """
        for i in range(self.number):
            i_rc = self.dna_to_binary(Seq.reverse_complement(
                self.binary_to_dna(i)))

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
        versa.

        :return: The forward and reverse complement counts.
        :rtype: tuple(list[float], list[float])
        """
        forward = []
        reverse = []

        for i in range(self.number):
            i_rc = self.dna_to_binary(Seq.reverse_complement(
                self.binary_to_dna(i)))

            if i <= i_rc:
                # Todo: Palindromes end up in both list, effectively increasing
                #   the combined total. Not sure what the alternative is, in
                #   general we cannot divide the counts over the two list (odd
                #   counts).
                #   It is clear that we can never really know what to do with
                #   palindrome counts. We just have to decide on a behaviour
                #   here. Perhaps we want to maintain the invariant that the
                #   combined total is the same as the original total? Then we
                #   could just divide by two and in case of an odd count put it
                #   in either one (always the first, random, or depending on
                #   the k-mer).
                forward.append(self.counts[i])
                reverse.append(self.counts[i_rc])
            #if
        #for

        return forward, reverse
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

        # Todo: store dtype='int64' in a module variable or something.
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

        :arg sequence:
        :type sequence: string

        :return: Binary representation of a DNA string.
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

        :arg number: Binary representation of a DNA string.
        :type number: int

        :returns A DNA string.
        :rtype: string
        """
        sequence = ""

        for i in range(self.length):
            sequence += self.__binary_to_nucleotide[number & 0x03]
            number >>= 2
        #while

        return sequence[::-1]
    #binary_to_dna

    def complement(self, number):
        """
        Calculate the complement of an integer (note, not the reverse
        complement).

        :arg number: An integer.
        :type number: int

        :return: The complement of {number}.
        :rtype: int
        """
        temp_number = number
        result = 0x00

        for i in range(self.length):
            result <<= 2
            result |= ~temp_number & 0x03
            temp_number <<= 2
        #for

        return result
    #complement

    def ratios_matrix(self):
        """
        Calculate all relative frequencies of k-mers. If a division by 0
        occurs, the frequency will be set to -1.0.

        :return: A matrix with relative frequencies.
        :rtype: float[][]
        """
        # Initialise the matrix.
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

    def freq_diff_matrix(self):
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

    def print_ratios(self, ratios):
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
#kMer
