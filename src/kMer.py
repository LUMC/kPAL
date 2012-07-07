#!/usr/bin/python

"""
Make a k-mer profile from a fastq or a fasta file.


"""

import sys
import argparse
from Bio import SeqIO

class kMer() :
    """
    Handle the counting of k-mers in a fastq file.
    """

    def __init__(self, length, inputFormat="fastq") :
        """
        Initialise the class instance.

        @arg length: Length of the k-mers.
        @type length: integer
        @arg inputFormat: Input format, either fastq or fasta.
        @type inputFormat: str
        """
        self.totalKMers = 0
        self.nonZeroKMers = 0
        self.inputFormat = inputFormat

        if length :
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
        for i in self.__nucleotideToBinary :
            self.__binaryToNucleotide[self.__nucleotideToBinary[i]] = i
    #__init__

    def __initialise(self, length) :
        """
        Initialise the class insance on demand (used by either __init__() or
        loadKMerCounts().

        @arg length: Length of the k-mers.
        @type length: integer
        """
        self.kMerLength = length
        self.numberOfKMers = 4 ** self.kMerLength
        self.kMerCount = self.numberOfKMers * [0]
        self.__bitMask = self.numberOfKMers - 0x01
    #__initialise

    def __addKMer(self, binaryRepresentation) :
        """
        Add a k-mer and keep track of the nonZeroKMers and totalKMers counters.

        @arg binaryRepresentation: Binary representation of a k-mer.
        @type binaryRepresentation: integer
        """
        if not self.kMerCount[binaryRepresentation] :
            self.nonZeroKMers += 1
        self.kMerCount[binaryRepresentation] += 1
        self.totalKMers += 1
    #__addKMer

    def __scanLine(self, sequence) :
        """
        Count all occurrences of  k-mers in a DNA string.

        @arg sequence: A DNA sequence from a fasta/fastq file.
        @type sequence: string
        """
        if len(sequence) > self.kMerLength :
            binaryRepresentation = 0x00

            # Calculate the binary representation of a k-mer.
            for i in sequence[:self.kMerLength] :
                binaryRepresentation = ((binaryRepresentation << 2) |
                    self.__nucleotideToBinary[i])
            self.__addKMer(binaryRepresentation)

            # Calculate the binary representation of the next k-mer.
            for i in sequence[self.kMerLength:] :
                binaryRepresentation = ((binaryRepresentation << 2) |
                    self.__nucleotideToBinary[i]) & self.__bitMask
                self.__addKMer(binaryRepresentation)
            #for
        #if
    #__scanLine

    def scanFastq(self, handle) :
        """
        Read a fasta/fastq file and count all k-mers in each line.

        @arg handle: An open file handle to a fastq file.
        @type handle: stream
        """
        for record in SeqIO.parse(handle, self.inputFormat) :
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

        for i in range(self.kMerLength) :
            sequence += self.__binaryToNucleotide[number & 0x03]
            number >>= 2
        #while

        return sequence[::-1]
    #binaryToDNA

    def saveKMerCounts(self, handle) :
        """
        Save the k-mer table in a file.

        @arg handle: Writable open handle to a file.
        @type handle: stream
        """
        handle.write("%i\n" % self.kMerLength)
        handle.write("%i\n" % self.totalKMers)
        handle.write("%i\n" % self.nonZeroKMers)
        for i in self.kMerCount :
            handle.write("%i\n" % i)
    #saveKMerCounts

    def loadKMerCounts(self, handle) :
        """
        Load the k-mer table from a file.

        @arg handle: Open handle to a file.
        @type handle: stream
        """
        self.__initialise(int(handle.readline()[:-1]))
        self.totalKMers = (int(handle.readline()[:-1]))
        self.nonZeroKMers = (int(handle.readline()[:-1]))

        offset = 0
        line = handle.readline()
        while line :
            self.kMerCount[offset] = int(line[:-1])
            line = handle.readline()
            offset += 1
        #while
    #loadKMerCounts

    def mergeKMerCounts(self, kMerInstance) :
        """
        Add the counts of a (compatible) kMer instance to this one.

        @arg kMerInstance: An other kMer instance.
        @type kMerInstance: kMer
        """
        self.totalKMers += kMerInstance.totalKMers
        self.nonZeroKMers = 0

        for i in range(self.numberOfKMers) :
            self.kMerCount[i] += kMerInstance.kMerCount[i]
            if self.kMerCount[i] :
                self.nonZeroKMers += 1
        #for
    #mergeKMerCounts

    def calculateRatios(self) :
        """
        Calculate all relative frequencies of k-mers. If a division by 0
        occurs, the frequency will be set to -1.0.

        @returns: A matrix with relative frequencies.
        @rtype: float[][]
        """
        # Initialise the matrix.
        ratios = []
        for i in range(self.numberOfKMers) :
            ratios.append(self.numberOfKMers * [0.0])

        # Fill the matrix.
        for i in range(self.numberOfKMers) :
            for j in range(self.numberOfKMers) :
                if self.kMerCount[j] :
                    ratios[i][j] = (float(self.kMerCount[i]) /
                        self.kMerCount[j]) / self.totalKMers
                else :
                    ratios[i][j] = -1.0
            #for

        return ratios
    #calculateRatios

    def calculateDiff(self) :
        """
        Calculate all frequency differences of k-mers.

        @returns: A matrix with frequency differences.
        @rtype: float[][]
        """
        ratios = []
        for i in range(self.numberOfKMers) :
            ratios.append(self.numberOfKMers * [0])

        # Fill the matrix.
        for i in range(self.numberOfKMers) :
            for j in range(self.numberOfKMers) :
                if self.kMerCount[j] :
                    ratios[i][j] = float(abs(self.kMerCount[i]  -
                        self.kMerCount[j])) / self.totalKMers

        return ratios
    #calculateDiff

    def printCounts(self) :
        """
        Print the k-mer counts.
        """
        for i in range(self.numberOfKMers) :
            print self.binaryToDNA(i), self.kMerCount[i]
    #printCounts

    def printRatios(self, ratios) :
        """
        Print a ratios matrix.

        @arg ratios: A matrix with relative frequencies.
        @type ratios: float[][]
        """
        # The header.
        print (self.kMerLength + 1) * ' ',
        for i in range(self.numberOfKMers) :
            print "%s%s" % (self.binaryToDNA(i), 2 * ' '),
        print

        # The matrix.
        for i in range(self.numberOfKMers) :
            print "%s" % self.binaryToDNA(i),
            for j in range(self.numberOfKMers) :
                print ("%%.%if" % self.kMerLength) % ratios[i][j],
            print
        #for
    #printRatios
#kMer

def makeProfile(input, output, kMerSize, inputFormat) :
    """
    """
    kMerInstance = kMer(kMerSize, inputFormat)
    kMerInstance.scanFastq(input)
    kMerInstance.saveKMerCounts(output)
#makeProfile

def main() :
    """
    Main entry point.
    """
    usage = __doc__.split("\n\n\n")
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=usage[0],
        epilog=usage[1])

    parser.add_argument("-i", dest="input", type=argparse.FileType('r'),
        required=True, help="the input file in fasta/fastq format")
    parser.add_argument("-o", dest="output", type=argparse.FileType('w'),
        required=True, help="the output file name")
    parser.add_argument("-a", dest="inputFormat", default="fastq",
        action="store_const", const="fasta",
        help="use fasta instead of fastq as input format")
    parser.add_argument("-k", dest="kMerSize", type=int, required=True,
        help="size of the k-mers")

    arguments = parser.parse_args()

    makeProfile(arguments.input, arguments.output, arguments.kMerSize,
        arguments.inputFormat)
#main

if __name__ == "__main__" :
    main()
