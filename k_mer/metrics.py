#!/usr/bin/python

"""
General library containing metrics and helper functions.
"""

import collections
import math

def median(vector):
    """
    Calculate the median of a vector.

    :arg vector: A vector.
    :type vector: list[float]

    :return: The median of the values in {vector}.
    :rtype: float
    """
    sorted_vector = sorted(vector)
    length = len(vector)

    if length % 2:
        position = length / 2
        return (sorted_vector[position] + sorted_vector[position + 1]) / 2.0
    #if
    return sorted_vector[length / 2]
#median

def distribution(vector):
    """
    Calculate the distribution of the values in a vector.

    :arg vector: A vector.
    :type vector: list[int]

    :return: A list of pairs (count, occurrences).
    :rtype: list[tuple(int, in)]
    """
    d = collections.defaultdict(int)

    for i in vector:
        d[i] += 1
    return map(lambda x: (x, d[x]), sorted(d))
#distribution

def vector_length(vector):
    """
    Calculate the Euclidean length of a vector.

    :arg vector: A vector.
    :type vector: list[int]

    :return: The length of {vector}.
    :rtype: float
    """
    return math.sqrt(sum(map(lambda x: x * x, vector)))
#vector_length

def stats(l):
    """
    Calculate the mean and median of l.

    :arg l: A list of values.
    :type l: list(float)

    :return: The mean and standard deviation of l.
    :rtype: tuple(float, float)
    """
    sum_l = 0
    sum_squared_l = 0
    n = 0

    for i in l:
        sum_l += i
        sum_squared_l += i * i
        n += 1
    #for
    
    mean = sum_l / float(n)
    return mean, math.sqrt((sum_squared_l / float(n)) - (mean * mean))
#stats

def get_scale(vector1, vector2):
    """
    Calculate scaling factors based upon total counts. One of the factors
    is always one (the other is either one or larger than one).

    :arg vector1: A vector.
    :type vector1: list[int]
    :arg vector2: A vector.
    :type vector2: list[int]

    :return: A tuple of scaling factors.
    :rtype: tuple(float)
    """
    scale1 = 1.0
    scale2 = 1.0

    # Calculate the scaling factors in such a way that no element is
    # between 0 and 1.
    vector1_total = sum(vector1)
    vector2_total = sum(vector2)

    if vector1_total < vector2_total:
        scale1 = float(vector2_total) / vector1_total
    else:
        scale2 = float(vector1_total) / vector2_total

    return scale1, scale2
#get_scale

def scale_down(scale1, scale2):
    """
    Normalise scaling factor between 0 and 1.

    :arg scale1: Scaling factor.
    :type scale1: float
    :arg scale2: Scaling factor.
    :type scale2: float

    :return: Tuple of normalised scaling factors.
    :rtype: tuple(float, float)
    """
    factor = max(scale1, scale2)

    return scale1 / factor, scale2 / factor
#scale_down

def scale(vector, scale):
    """
    Scale a vector.

    :arg vector: A vector.
    :type vector: list[float]
    :arg scale: Scaling factor.
    :type scale: float

    :return: A vector scaled by factor {scale}.
    :rtype: list[float]
    """
    return map(lambda x: scale * x, vector)
#scale

def positive(vector1, vector2):
    """
    Set all zero positions in {vector2} to zero in {vector1}.

    :arg vector1: A vector.
    :type vector1: list[int]
    :arg vector2: A vector.
    :type vector2: list[int]

    :return: {vector1} with all zero positions in {vector2} set to zero.
    :rtype: list[int]
    """
    return map(lambda x: x[1] and x[0], zip(vector1, vector2))
#positive

def multiset(vector1, vector2, pairwise):
    """
    Calculate the multiset distance between two vectors.

    :arg vector1: A vector.
    :type vector1: list[float]
    :arg vector2: A vector.
    :type vector2: list[float]
    :arg pairwise: A pairwise distance function.
    :type pairwise: function

    :return: The multiset distance between {vector1} and {vector2}.
    :rtype: float
    """
    c = 0.0
    d = 1

    # Calculate the counter and the denominator of the distance function.
    for i in range(len(vector1)):
        if vector1[i] or vector2[i]:
            c += pairwise(vector1[i], vector2[i])
            d += 1
        #if
    #for

    return c / d
#multiset

def euclidean(vector1, vector2):
    """
    Calculate the Euclidean distance between two vectors.

    :arg vector1: A vector.
    :type vector1: list[float]
    :arg vector2: A vector.
    :type vector2: list[float]

    :return: The Euclidean distance between {vector1} and {vector2}.
    :rtype: float
    """
    sum_of_squares = 1

    # Calculate the counter and the denominator of the distance function.
    for i in range(len(vector1)):
        sum_of_squares += (vector1[i] - vector2[i]) ** 2

    return math.sqrt(sum_of_squares)
#euclidean

def cosine_similarity(vector1, vector2):
    """
    Calculate the Cosine similarity between two vectors.

    :arg vector1: A vector.
    :type vector1: list[float]
    :arg vector2: A vector.
    :type vector2: list[float]

    :return: The Cosine similarity between {vector1} and {vector2}.
    :rtype: float
    """
    return (sum(map(lambda x: x[0] * x[1], zip(vector1, vector2))) /
        (vector_length(vector1) * vector_length(vector2)))
#cosine_similarity

vector_distance = {
    "default": None,
    "euclidean": euclidean,
    "cosine": cosine_similarity
}
""" Dictionary of vector distance functions. """

pairwise = {
    "diff-prod": lambda x, y: abs(x - y) / float((x + 1) * (y + 1)),
    "diff-sum": lambda x, y: abs(x - y) / float(x + y + 1)
}
""" Dictionary of pairwise distance functions. """

summary = {
    "min": min,
    "average": lambda x: sum(x) / float(len(x)),
    "median": median
}
""" Dictionary of summary functions. """

mergers = {
    "sum": lambda x, y: x + y,
    "xor": lambda x, y: (x + y) * (not (x and y)),
    "int": lambda x, y: x * bool(y),
    "nint": lambda x, y: x * (not bool(y))
}
""" Merge functions. """
