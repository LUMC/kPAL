#!/usr/bin/python

"""
Toolbox for k-mer profiles.


"""

import argparse
import kLib
import kDiffLib

def makeProfile(input, output, size, inputFormat):
    """
    Make a k-mer profile from a fastq or a fasta file.

    @arg input: Open readable stream to a fasta or fastq file.
    @type input: stream
    @arg output: Open writeable stream.
    @type output: stream
    @arg size: k-mer size.
    @type size: int
    @arg inputFormat: Input format, either "fastq" or "fasta".
    @type inputFormat: str
    """
    profile = kLib.kMer(size, inputFormat)
    profile.scanFastq(input)
    profile.save(output)
#makeProfile

def mergeProfiles(input1, input2, output):
    """
    Merge two k-mer profiles.

    @arg input1: Open readable stream to a k-mer file.
    @type input1: stream
    @arg input2: Open readable stream to a k-mer file.
    @type input2: stream
    @arg input2: Open writeable stream to a k-mer file.
    @type input2: stream
    """
    kMerIn = kLib.kMer(0)
    kMerOut = kLib.kMer(0)

    kMerIn.load(input1)
    kMerOut.load(input2)

    if kMerIn.length != kMerOut.length:
        raise ValueError(kMer.lengthError)

    kMerOut.merge(kMerIn)
    kMerOut.save(output)
#mergeProfiles

def profileDiff(input1, input2, algorithm, precision, down):
    """
    Calculate the difference between two k-mer profiles.

    @arg input1: Open readable stream to a k-mer file.
    @type input1: stream
    @arg input2: Open readable stream to a k-mer file.
    @type input2: stream
    @arg algorithm: Selector for the distance algorithm.
    @type algorithm: int
    @arg precision: Number of digits in the output.
    @type precision: int
    @arg down: Selector for the scaling function.
    @type down: bool
    """
    kDiff = kDiffLib.kMerDiff(algorithm, down=down)
    if not kDiff.distance:
        print kDiffLib.kMerDiff.algorithmError
        parser.print_usage()
        return
    #if

    profile1 = kLib.kMer(0)
    profile2 = kLib.kMer(0)

    profile1.load(input1)
    profile2.load(input2)

    if profile1.length != profile2.length:
        raise ValueError(kLib.lengthError)

    return ("%%.%if" % precision) % kDiff.distance(profile1,
        profile2)
#profileDiff

def diffMatrix(inputs, output, algorithm, precision, down):
    """
    Make a distance matrix any number of k-mer profiles.

    @arg inputs: List of open readable streams to k-mer files.
    @type inputs: list[stream]
    @arg output: Open handle to a writable file.
    @type output: stream
    @arg algorithm: Selector for the distance algorithm.
    @type algorithm: int
    @arg precision: Number of digits in the output.
    @type precision: int
    @arg down: Selector for the scaling function.
    @type down: bool
    """
    if len(inputs) < 2:
        raise ValueError("You must give at least two input files.")

    kDiff = kDiffLib.kMerDiff(algorithm, down=down)
    if not kDiff.distance:
        raise ValueError(kDiffLib.kMerDiff.algorithmError)

    counts = []
    for i in inputs:
        counts.append(kLib.kMer(0))
        counts[-1].load(i)
        if counts[0].length != counts[-1].length:
            raise ValueError(kLib.lengthError)
    #for

    kDiffLib.makeDistanceMatrix(counts, output, precision, kDiff)
#diffMatrix

def smoothProfiles(input1, input2, output1, output2, function, threshold):
    """
    Dynamically smooth two k-mer profiles.

    @arg input1: Open readable stream to a k-mer file.
    @type input1: stream
    @arg input2: Open readable stream to a k-mer file.
    @type input2: stream
    @arg output1: Open writeable stream to a k-mer file.
    @type output1: stream
    @arg output2: Open writeable stream to a k-mer file.
    @type output2: stream
    @arg function: Function to assess the quality of a sub-profile.
    @type function: function
    @arg threshold: Threshold under which a sub-profile is collapsed.
    @type threshold: int
    """
    profile1 = kLib.kMer(0)
    profile2 = kLib.kMer(0)

    profile1.load(input1)
    profile1.load(input2)

    if profile1.length != profile2.length:
        raise ValueError(kLib.lengthError)

    dynamicSmooth(profile1, profile2, 0, profile1.number, function, threshold)

    profile1.save(output1)
    profile2.save(output2)
#smoothProfiles

#
# Wrapper functions for the command line interface.
#

def index(arguments):
    makeProfile(arguments.input, arguments.output, arguments.size,
        arguments.inputFormat)

def merge(arguments):
    mergeProfiles(arguments.input[0], arguments.input[1], arguments.output)

def diff(arguments):
    print profileDiff(arguments.input[0], arguments.input[1],
        arguments.algorithm, arguments.precision, arguments.down)

def matrix(arguments):
    diffMatrix(arguments.inputs, arguments.output, arguments.algorithm,
        arguments.precision, arguments.down)

def smooth(arguments):
    smoothProfiles(arguments.input[0], arguments.input[1], arguments.output[0],
        arguments.output[1], arguments.function, arguments.threshold)

def main():
    """
    Main entry point.
    """
    usage = __doc__.split("\n\n\n")
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=usage[0], epilog=usage[1])
    subparsers = parser.add_subparsers()

    parser_index = subparsers.add_parser("index",
        description=makeProfile.__doc__.split("\n\n")[0])
    parser_index.add_argument("-i", dest="input", type=argparse.FileType('r'),
        required=True, help="input file in fasta/fastq format")
    parser_index.add_argument("-o", dest="output", type=argparse.FileType('w'),
        required=True, help="output file name")
    parser_index.add_argument("-a", dest="inputFormat", default="fastq",
        action="store_const", const="fasta",
        help="use fasta instead of fastq as input format")
    parser_index.add_argument("-k", dest="size", type=int, required=True,
        help="k-mer size")
    parser_index.set_defaults(func=index)

    parser_merge = subparsers.add_parser("merge",
        description=mergeProfiles.__doc__.split("\n\n")[0])
    parser_merge.add_argument("-i", dest="input", type=argparse.FileType('r'),
        nargs=2, help="input file names")
    parser_merge.add_argument("-o", dest="output", type=argparse.FileType('w'),
        required=True, help="output file name")
    parser_merge.set_defaults(func=merge)

    parser_diff = subparsers.add_parser("diff",
        description=profileDiff.__doc__.split("\n\n")[0])
    parser_diff.add_argument("-i", dest="input", type=argparse.FileType('r'),
        required=True, nargs=2, help="input file names")
    parser_diff.add_argument("-p", dest="precision", type=int, default=3, 
        help="number of decimals")
    parser_diff.add_argument("-a", dest="algorithm", type=int, default=0,
        help="distance algorithm")
    parser_diff.add_argument("-d", dest="down", default=False,
        action="store_true", help="scale down")
    parser_diff.set_defaults(func=diff)

    parser_matrix = subparsers.add_parser("matrix",
        description=diffMatrix.__doc__.split("\n\n")[0])
    parser_matrix.add_argument("-i", dest="inputs", nargs='+',
        type=argparse.FileType('r'), help="input file names")
    parser_matrix.add_argument("-o", dest="output", required=True,
        type=argparse.FileType('w'), help="output file name")
    parser_matrix.add_argument("-p", dest="precision", type=int,
        default=3, help="number of decimals")
    parser_matrix.add_argument("-a", dest="algorithm", type=int, default=0,
        help="distance algorithm")
    parser_matrix.add_argument("-d", dest="down", default=False,
        action="store_true", help="scale down")
    parser_matrix.set_defaults(func=matrix)

    parser_smooth = subparsers.add_parser("smooth",
        description=smoothProfiles.__doc__.split("\n\n")[0])
    parser_smooth.add_argument("-i", dest="input", type=argparse.FileType('r'),
        required=True, nargs=2, help="input file names")
    parser_smooth.add_argument("-o", dest="output", required=True, nargs=2,
        type=argparse.FileType('w'), help="output file names")
    parser_smooth.set_defaults(func=smooth)

    arguments = parser.parse_args()

    arguments.func(arguments)
#main

if __name__ == "__main__":
    main()
