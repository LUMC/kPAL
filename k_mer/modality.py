#!/usr/bin/env python

import argparse
import math
import pylab

from scipy import interpolate

def modality(handle, cutoff, resolution):
    """
    """
    x, y = zip(*map(lambda x: map(int, x.split()),
        handle.readlines()[:cutoff]))

    # Ad-hoc weighing function.
    w = map(lambda x: 1.0 / (10 * x + 1), x)
    # Get an approximation of the distribution.
    f = interpolate.UnivariateSpline(x, y, w=w)

    xs = pylab.linspace(0, cutoff, resolution)

    # Plot the original.
    pylab.plot(x, y)

    # Plot the interpolated function.
    ys = f(xs)
    pylab.plot(xs, ys)

    # Plot the tops.
    ys = f(xs, 1)
    g = interpolate.InterpolatedUnivariateSpline(xs, ys)
    pylab.plot(g.roots(), f(g.roots()), 'o')

    # Plot the bending points.
    ys = f(xs, 2)
    g = interpolate.InterpolatedUnivariateSpline(xs, ys)
    pylab.plot(g.roots(), f(g.roots()), 'o')

    pylab.legend(('original', 'interpolated', '1st derivative',
        '2nd derivative'), loc='best')
    pylab.xscale("log")
    pylab.yscale("log")
    pylab.show()
#modality

def main():
    """
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("INPUT", type=argparse.FileType('r'),
        help="distribution file")
    parser.add_argument("-c", dest="cutoff", type=int, default=150,
        help="cutoff value")
    parser.add_argument("-r", dest="resolution", type=int, default=10000,
        help="resolution")

    args = parser.parse_args()

    modality(args.INPUT, args.cutoff, args.resolution)
#main

if __name__ == "__main__":
    main()
