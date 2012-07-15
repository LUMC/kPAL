#!/usr/bin/python

"""
k-mer base library.
"""

import sys
import argparse
from Bio import SeqIO

class kMer():
    """
    Handle the counting of k-mers in a fastq file.
    """
    lengthError = "k-mer lengths of the files differ."

    def __init__(self, length, inputFormat="fastq"):
        """
        Initialise the class instance.

        @arg length: Length of the k-mers.
        @type length: integer
        @arg inputFormat: Input format, either fastq or fasta.
        @type inputFormat: str
        """
        self.total = 0
        self.nonZero = 0
        self.inputFormat = inputFormat

        if length:
            self.__initialise(length)

        # Conversion table form nucleotides to binary.
        self.__nucleotideToBinary = {
            'A' : 0x00,
            'C' : 0x01,
            'G' : 0x02,
            'T' : 0x03
        }

        # Build the reverse table of __nucleotideToBinary.
        self.__binaryToNucleotide = {}
        for i in self.__nucleotideToBinary:
            self.__binaryToNucleotide[self.__nucleotideToBinary[i]] = i
    #__init__

    def __initialise(self, length):
        """
        Initialise the class insance on demand (used by either __init__() or
        load().

        @arg length: Length of the k-mers.
        @type length: integer
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
        @type binaryRepresentation: integer
        """
        if not self.count[binaryRepresentation]:
            self.nonZero += 1
        self.count[binaryRepresentation] += 1
        self.total += 1
    #__add

    def __scanLine(self, sequence):
        """
        Count all occurrences of  k-mers in a DNA string.

        @arg sequence: A DNA sequence from a fasta/fastq file.
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

    def scanFastq(self, handle):
        """
        Read a fasta/fastq file and count all k-mers in each line.

        @arg handle: An open file handle to a fastq file.
        @type handle: stream
        """
        for record in SeqIO.parse(handle, self.inputFormat):
            for sequence in str(record.seq).split('N'):
                self.__scanLine(sequence)
    #scanFastq

    def binaryToDNA(self, number):
        """
        Convert an integer to a DNA string.

        @arg number: Binary representation of a DNA string.
        @type number: integer

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

    def load(self, handle):
        """
        Load the k-mer table from a file.

        @arg handle: Open handle to a file.
        @type handle: stream
        """
        self.__initialise(int(handle.readline()[:-1]))
        self.total = (int(handle.readline()[:-1]))
        self.nonZero = (int(handle.readline()[:-1]))

        offset = 0
        line = handle.readline()
        while line:
            self.count[offset] = int(line[:-1])
            line = handle.readline()
            offset += 1
        #while
    #load

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
