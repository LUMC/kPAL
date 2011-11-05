#!/usr/bin/python

"""
@requires: sys
@requires: Bio.SeqIO
"""

import sys
from Bio import SeqIO

class kMer() :
    """
    """

    def __init__(self, length) :
        """
        @arg length: Length of the k-mers.
        @type length: integer
        """

        self.__kMerLength = length
        self.__numberOfKMers = 4 ** self.__kMerLength
        self.__bitMask = self.__numberOfKMers - 1
        self.__kMerCount = self.__numberOfKMers * [0]

        # Conversion table form nucleotides to binary.
        self.__nucleotideToBinary = {
            'A' : 0,
            'C' : 1,
            'G' : 2,
            'T' : 3
        }

        # Build the reverse table of __nucleotideToBinary.
        self.__binaryToNucleotide = {}
        for i in self.__nucleotideToBinary :
            self.__binaryToNucleotide[self.__nucleotideToBinary[i]] = i
    #__init__

    def __scanLine(self, sequence) :
        """
        Count all occurrences of  k-mers in a DNA string.

        @arg sequence: A DNA sequence from a fastq file.
        @type sequence: string
        """

        if len(sequence) > self.__kMerLength :
            binaryRepresentation = 0

            # Calculate the binary representation of a k-mer.
            for i in sequence[:self.__kMerLength] :
                binaryRepresentation = ((binaryRepresentation << 0x02) |
                    self.__nucleotideToBinary[i])
            self.__kMerCount[binaryRepresentation] += 1

            # Calculate the binary representation of the next k-mer.
            for i in sequence[self.__kMerLength:] :
                binaryRepresentation = ((binaryRepresentation << 0x02) |
                    self.__nucleotideToBinary[i]) & self.__bitMask
                self.__kMerCount[binaryRepresentation] += 1
            #for
        #if
    #__scanLine

    def scanFastq(self, handle) :
        """
        Read a fastq file and count all k-mers in each line.

        @arg handle: An open file handle to a fastq file.
        @type handle: stream
        """

        for record in SeqIO.parse(handle, "fastq") :
            for sequence in str(record.seq).split('N') :
                self.__scanLine(sequence)
    #scanFastq

    def binaryToDNA(self, number) :
        """
        Convert an integer to a DNA string.

        @arg number: Binary representation of a DNA string.
        @type number: integer

        @returns: A DNA string.
        @rtype: string
        """

        sequence = ""

        for i in range(self.__kMerLength) :
            sequence += self.__binaryToNucleotide[number & 0x03]
            number >>= 2
        #while

        return sequence[::-1]
    #binaryToDNA

    def calculateRatios(self) :
        """
        Calculate all relative frequencies of k-mers. If a division by 0
        occurs, the frequency will be set to -1.0.

        @returns: A matrix with relative frequencies.
        @rtype: float[][]
        """

        # Initialise the matrix.
        ratios = []
        for i in range(self.__numberOfKMers) :
            ratios.append(self.__numberOfKMers * [0.0])

        # Fill the matrix.
        for i in range(self.__numberOfKMers) :
            for j in range(self.__numberOfKMers) :
                if self.__kMerCount[j] :
                    ratios[i][j] = \
                        float(self.__kMerCount[i]) / self.__kMerCount[j]
                else :
                    ratios[i][j] = -1.0
            #for

        return ratios
    #calculateRatios

    def printCounts(self) :
        """
        Print the k-mer counts.
        """

        for i in range(self.__numberOfKMers) :
            print self.binaryToDNA(i), self.__kMerCount[i]
    #printCounts

    def printRatios(self, ratios) :
        """
        Print a ratios matrix.
        """

        # The header.
        print (self.__kMerLength + 1) * ' ',
        for i in range(self.__numberOfKMers) :
            print "%s%s" % (self.binaryToDNA(i), 2 * ' '),
        print

        # The matrix.
        for i in range(self.__numberOfKMers) :
            print "%s" % self.binaryToDNA(i),
            for j in range(self.__numberOfKMers) :
                print ("%%.%if" % self.__kMerLength) % ratios[i][j],
            print
        #for
    #printRatios
#kMer

def main() :
    """
    """

    kMerInstance = kMer(2)

    inputHandle = open(sys.argv[1], "r")
    kMerInstance.scanFastq(inputHandle)
    inputHandle.close()

    kMerInstance.printCounts()

    ratios = kMerInstance.calculateRatios()
    kMerInstance.printRatios(ratios)
#main

if __name__ == "__main__" :
    main()
