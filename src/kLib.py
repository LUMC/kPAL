#!/usr/bin/python

"""
k-mer base library.
"""

import sys
import argparse
import re
import random
from os import path
from Bio import Seq
from Bio import SeqIO

class kMer():
    """
    Handle the counting of k-mers in a fasta file.
    """
    __nucleotideToBinary = {
        'A': 0x00, 'a': 0x00,
        'C': 0x01, 'c': 0x01,
        'G': 0x02, 'g': 0x02,
        'T': 0x03, 't': 0x03
    }
    """ Conversion table form nucleotide to binary. """
    __binaryToNucleotide = {
        0x00: 'A',
        0x01: 'C',
        0x02: 'G',
        0x03: 'T'
    }
    """ Conversion table form binary to nucleotide. """

    def __init__(self):
        """
        Initialise the class instance.
        """
        self.total = 0
        self.nonZero = 0
        self.name = ""
    #__init__

    def __initialise(self, length):
        """
        Initialise the rest of the class instance variables.

        @arg length: Length of the k-mers.
        @type length: int
        """
        self.length = length
        self.number = 4 ** self.length
        self.count = self.number * [0]
        self.__bitMask = self.number - 0x01
    #__initialise

    def __add(self, binaryRepresentation):
        """
        Add a k-mer and keep track of the nonZero and total counters.

        @arg binaryRepresentation: Binary representation of a k-mer.
        @type binaryRepresentation: int
        """
        if not self.count[binaryRepresentation]:
            self.nonZero += 1
        self.count[binaryRepresentation] += 1
        self.total += 1
    #__add

    def __scanLine(self, sequence):
        """
        Count all occurrences of  k-mers in a DNA string.

        @arg sequence: A DNA sequence from a fasta file.
        @type sequence: string
        """
        if len(sequence) > self.length:
            binaryRepresentation = 0x00

            # Calculate the binary representation of a k-mer.
            for i in sequence[:self.length]:
                binaryRepresentation = ((binaryRepresentation << 2) |
                    self.__nucleotideToBinary[i])
            self.__add(binaryRepresentation)

            # Calculate the binary representation of the next k-mer.
            for i in sequence[self.length:]:
                binaryRepresentation = ((binaryRepresentation << 2) |
                    self.__nucleotideToBinary[i]) & self.__bitMask
                self.__add(binaryRepresentation)
            #for
        #if
    #__scanLine

    def analyse(self, handle, length):
        """
        Read a fasta file and count all k-mers in each line.

        @arg handle: An open file handle to a fasta file.
        @type handle: stream
        @arg length: Length of the k-mers.
        @type length: int
        """
        self.__initialise(length)
        S = re.compile("[^%s]" % ''.join(self.__nucleotideToBinary.keys()))

        for record in SeqIO.parse(handle, "fasta"):
            for sequence in S.split(str(record.seq)):
                self.__scanLine(sequence)
    #analyse

    def load(self, handle):
        """
        Load the k-mer table from a file.

        @arg handle: Open handle to a file.
        @type handle: stream
        """
        self.__initialise(int(handle.readline()[:-1]))
        self.total = (int(handle.readline()[:-1]))
        self.nonZero = (int(handle.readline()[:-1]))
        self.name = path.basename(handle.name)

        offset = 0
        line = handle.readline()
        while line:
            self.count[offset] = int(line[:-1])
            line = handle.readline()
            offset += 1
        #while
    #load

    def save(self, handle):
        """
        Save the k-mer table in a file.

        @arg handle: Writable open handle to a file.
        @type handle: stream
        """
        handle.write("%i\n" % self.length)
        handle.write("%i\n" % self.total)
        handle.write("%i\n" % self.nonZero)
        for i in self.count:
            handle.write("%i\n" % i)
    #save

    def merge(self, profile):
        """
        Add the counts of a (compatible) k-mer profile to this one.

        @arg profile: An other k-mer profile.
        @type profile: object(kMer)
        """
        self.total += profile.total
        self.nonZero = 0

        for i in range(self.number):
            self.count[i] += profile.count[i]
            if self.count[i]:
                self.nonZero += 1
        #for
    #merge

    def balance(self):
        """
        Add the counts of the reverse complement of a k-mer to the k-mer and
        vice versa.
        """
        for i in range(self.number):
            i_rc = self.DNAToBinary(Seq.reverse_complement(
                self.binaryToDNA(i)))

            if i <= i_rc:
                temp = self.count[i]
                self.count[i] += self.count[i_rc]
                self.count[i_rc] += temp
            #if
        #for
    #balance

    def split(self):
        """
        Split the profile into two lists, every position in the first list has
        its reverse complement in the same position in the second list and vice
        versa.

        @returns: The forward and reverse complement counts.
        @rtype: tuple(list[float], list[float])
        """
        forward = []
        reverse = []

        for i in range(self.number):
            i_rc = self.DNAToBinary(Seq.reverse_complement(
                self.binaryToDNA(i)))

            if i <= i_rc:
                forward.append(self.count[i])
                reverse.append(self.count[i_rc])
            #if
        #for

        return forward, reverse
    #split

    def shrink(self, factor):
        """
        Shrink the profile, effectively reducing the value of k.

        Note that this operation may give slightly different values than
        indexing on a lower k directly.

        @arg factor: Shrinking factor.
        @type factor: int
        """
        if self.length < factor:
            raise ValueError(
                "Reduction factor should be smaller than k-mer size.")

        mergeSize = 4 ** factor
        self.nonZero = 0
        newCount = []
        for i in range(0, self.number, mergeSize):
            sub = sum(map(lambda x: self.count[x], range(i, i + mergeSize)))

            if sub:
                self.nonZero += 1
            newCount.append(sub)
        #for
        self.count = newCount
        self.length -= factor
        self.number -= mergeSize
        self.__bitMask >>= 2 * factor
    #shrink

    def shuffle(self):
        """
        Randomise the profile.
        """
        random.shuffle(self.count)
    #shuffle

    def DNAToBinary(self, sequence):
        """
        Convert a string of DNA to an integer.

        @arg sequence:
        @type sequence: string

        @returns: Binary representation of a DNA string.
        @rtype: int
        """
        result = 0x00

        for i in sequence:
            result <<= 2
            result |= self.__nucleotideToBinary[i]
        #for

        return result
    #DNAToBinary

    def binaryToDNA(self, number):
        """
        Convert an integer to a DNA string.

        @arg number: Binary representation of a DNA string.
        @type number: int

        @returns: A DNA string.
        @rtype: string
        """
        sequence = ""

        for i in range(self.length):
            sequence += self.__binaryToNucleotide[number & 0x03]
            number >>= 2
        #while

        return sequence[::-1]
    #binaryToDNA

    def complement(self, number):
        """
        Calculate the complement of an integer (note, not the reverse
        complement).

        @arg number: An integer.
        @type number: int

        @returns: The complement of {number}.
        @rtype: int
        """
        myNumber = number
        result = 0x00

        for i in range(self.length):
            result <<= 2
            result |= ~myNumber & 0x03
            myNumber <<= 2
        #for

        return result
    #complement

    def calculateRatios(self):
        """
        Calculate all relative frequencies of k-mers. If a division by 0
        occurs, the frequency will be set to -1.0.

        @returns: A matrix with relative frequencies.
        @rtype: float[][]
        """
        # Initialise the matrix.
        ratios = []
        for i in range(self.number):
            ratios.append(self.number * [0.0])

        # Fill the matrix.
        for i in range(self.number):
            for j in range(self.number):
                if self.count[j]:
                    ratios[i][j] = (float(self.count[i]) /
                        self.count[j]) / self.total
                else:
                    ratios[i][j] = -1.0
            #for

        return ratios
    #calculateRatios

    def calculateDiff(self):
        """
        Calculate all frequency differences of k-mers.

        @returns: A matrix with frequency differences.
        @rtype: float[][]
        """
        ratios = []
        for i in range(self.number):
            ratios.append(self.number * [0])

        # Fill the matrix.
        for i in range(self.number):
            for j in range(self.number):
                if self.count[j]:
                    ratios[i][j] = float(abs(self.count[i] -
                        self.count[j])) / self.total

        return ratios
    #calculateDiff

    def printCounts(self):
        """
        Print the k-mer counts.
        """
        for i in range(self.number):
            print self.binaryToDNA(i), self.count[i]
    #printCounts

    def printRatios(self, ratios):
        """
        Print a ratios matrix.

        @arg ratios: A matrix with relative frequencies.
        @type ratios: float[][]
        """
        # The header.
        print (self.length + 1) * ' ',
        for i in range(self.number):
            print "%s%s" % (self.binaryToDNA(i), 2 * ' '),
        print

        # The matrix.
        for i in range(self.number):
            print "%s" % self.binaryToDNA(i),
            for j in range(self.number):
                print ("%%.%if" % self.length) % ratios[i][j],
            print
        #for
    #printRatios
#kMer
