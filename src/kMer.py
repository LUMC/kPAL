#!/usr/bin/python

"""
@requires: sys
@requires: Bio.SeqIO
"""

import sys
from Bio import SeqIO

kMerLength = 2
numberOfKMers = 4 ** kMerLength
bitMask = numberOfKMers - 1
bucket = numberOfKMers * [0]

nucleotideToBinary = {
    'A' : 0,
    'C' : 1,
    'G' : 2,
    'T' : 3
}

binaryToNucleotide = {}

def binaryToDNA(number) :
    """
    """

    sequence = ""

    for i in range(kMerLength) :
        sequence += binaryToNucleotide[number & 0x03]
        number >>= 2
    #while

    return sequence[::-1]
#binaryToDNA

def kCount(sequence) :
    """
    """

    for DNASeq in sequence.split('N') :
        if len(DNASeq) > kMerLength :
            binaryRepresentation = 0
            for i in DNASeq[:kMerLength] :
                binaryRepresentation = ((binaryRepresentation << 0x02) |
                    nucleotideToBinary[i])
            bucket[binaryRepresentation] += 1
            for i in DNASeq[kMerLength:] :
                binaryRepresentation = ((binaryRepresentation << 0x02) |
                    nucleotideToBinary[i]) & bitMask
                bucket[binaryRepresentation] += 1
            #for
        #if
    #for
#kCount

def calculateRatios() :
    """
    """

    ratios = []
    for i in range(numberOfKMers) :
        ratios.append(numberOfKMers * [0.0])

    for i in range(numberOfKMers) :
        for j in range(numberOfKMers) :
            if bucket[j] :
                ratios[i][j] = float(bucket[i]) / bucket[j]
            else :
                ratios[i][j] = -1.0
        #for

    return ratios
#calculateRatios

def printRatios(ratios) :
    """
    """

    print (kMerLength + 1) * ' ',
    for i in range(numberOfKMers) :
        print "%s%s" % (binaryToDNA(i), 2 * ' '),
    print
    for i in range(numberOfKMers) :
        print "%s" % binaryToDNA(i),
        for j in range(numberOfKMers) :
            print ("%%.%if" % kMerLength) % ratios[i][j],
        print
    #for
#printRatios

def main() :
    """
    """

    for i in nucleotideToBinary :
        binaryToNucleotide[nucleotideToBinary[i]] = i

    inputHandle = open(sys.argv[1], "r")

    for record in SeqIO.parse(inputHandle, "fastq") :
        kCount(str(record.seq))

    #for i in range(numberOfKMers) :
    #    print binaryToDNA(i), bucket[i]

    ratios = calculateRatios()
    printRatios(ratios)

    inputHandle.close()
#main

if __name__ == "__main__" :
    main()
