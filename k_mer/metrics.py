"""
General library containing metrics and helper functions.
"""


from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future import standard_library

import numpy as np

with standard_library.hooks():
    from collections import Counter


def distribution(vector):
    """
    Calculate the distribution of the values in a vector.

    :arg vector: A vector.
    :type vector: iterable of int

    :return: A list of (value, count) pairs.
    :rtype: list of (int, int)
    """
    # Todo: I'm not sure this should be in this module.
    return sorted(Counter(vector).items())
#distribution


def vector_length(vector):
    """
    Calculate the Euclidean length of a vector.

    :arg vector: A vector.
    :type vector: array_like, 1 dimension

    :return: The length of {vector}.
    :rtype: float
    """
    # Note: This seems faster than `numpy.linalg.norm`.
    return np.sqrt(np.dot(vector, vector))
#vector_length


def get_scale(vector1, vector2):
    """
    Calculate scaling factors based upon total counts. One of the factors
    is always one (the other is either one or larger than one).

    :arg vector1, vector2: Vectors.
    :type vector1, vector2: array_like, 1 dimension

    :return: A tuple of scaling factors.
    :rtype: (float, float)
    """
    scale1 = 1.0
    scale2 = 1.0

    # Calculate the scaling factors in such a way that no element is
    # between 0 and 1.
    vector1_total = np.sum(vector1)
    vector2_total = np.sum(vector2)

    if vector1_total < vector2_total:
        scale1 = vector2_total / vector1_total
    else:
        scale2 = vector1_total / vector2_total

    return scale1, scale2
#get_scale


def scale_down(scale1, scale2):
    """
    Normalise scaling factor between 0 and 1.

    :arg scale1, scale2: Scaling factors.
    :type scale1, scale2: float

    :return: Tuple of normalised scaling factors.
    :rtype: (float, float)
    """
    factor = max(scale1, scale2)

    return scale1 / factor, scale2 / factor
#scale_down


def positive(vector1, vector2):
    """
    Set all zero positions in {vector2} to zero in {vector1}.

    :arg vector1, vector2: Vectors.
    :type vector1, vector2: array_like, 1 dimension

    :return: {vector1} with all zero positions in {vector2} set to zero.
    :rtype: ndarray
    """
    return np.multiply(vector1, np.asanyarray(vector2, dtype=bool))
#positive


def multiset(vector1, vector2, pairwise):
    """
    Calculate the multiset distance between two vectors.

    :arg vector1, vector2: Vectors.
    :type vector1, vector2: array_like, 1 dimension
    :arg pairwise: A pairwise distance function.
    :type pairwise: function

    :return: The multiset distance between {vector1} and {vector2}.
    :rtype: float

    Note that `function` must be vectorized, i.e., it is called directly on
    NumPy arrays, instead of on their pairwise elements. If your function only
    works on individual elements, convert it to a NumPy ufunc first. For
    example:

        f = np.vectorize(f, otypes=['float'])
    """
    vector1 = np.asanyarray(vector1)
    vector2 = np.asanyarray(vector2)

    nonzero = np.where(np.logical_or(vector1, vector2))
    distances = pairwise(vector1[nonzero], vector2[nonzero])
    return distances.sum() / (len(distances) + 1)
#multiset


def euclidean(vector1, vector2):
    """
    Calculate the Euclidean distance between two vectors.

    :arg vector1, vector2: Vectors.
    :type vector1, vector2: array_like, 1 dimension

    :return: The Euclidean distance between {vector1} and {vector2}.
    :rtype: float
    """
    return vector_length(np.subtract(vector2, vector1))
#euclidean


def cosine_similarity(vector1, vector2):
    """
    Calculate the Cosine similarity between two vectors.

    :arg vector1, vector2: Vectors.
    :type vector1, vector2: array_like, 1 dimension

    :return: The Cosine similarity between {vector1} and {vector2}.
    :rtype: float
    """
    return np.dot(vector1, vector2) / (vector_length(vector1) *
                                       vector_length(vector2))
#cosine_similarity


vector_distance = {
    "default": None,
    "euclidean": euclidean,
    "cosine": cosine_similarity
}
""" Dictionary of vector distance functions. """


# Note: These definitions are safe to call on arrays.
pairwise = {
    "diff-prod": lambda x, y: abs(x - y) / ((x + 1) * (y + 1)),
    "diff-sum": lambda x, y: abs(x - y) / (x + y + 1)
}
""" Dictionary of pairwise distance functions. """


summary = {
    "min": np.min,
    "average": np.mean,
    "median": np.median
}
""" Dictionary of summary functions. """


# Note: These definitions are safe to call on arrays.
mergers = {
    "sum": lambda x, y: x + y,
    "xor": lambda x, y: (x + y) * np.logical_xor(x, y),
    "int": lambda x, y: x * np.asanyarray(y, dtype=bool),
    "nint": lambda x, y: x * np.logical_not(y)
}
""" Merge functions. """
