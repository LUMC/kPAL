"""
General library containing metrics and helper functions.

.. moduleauthor:: Leiden University Medical Center <humgen@lumc.nl>
.. moduleauthor:: Jeroen F.J. Laros <j.f.j.laros@lumc.nl>
.. moduleauthor:: Martijn Vermaat <m.vermaat.hg@lumc.nl>

.. Licensed under the MIT license, see the LICENSE file.
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
    :type vector: iterable(int)

    :return: A list of `(value, count)` pairs.
    :rtype: list(int, int)
    """
    # Todo: I'm not sure this should be in this module.
    return sorted(Counter(vector).items())


def vector_length(vector):
    """
    Calculate the Euclidean length of a vector.

    :arg array_like vector: A vector.

    :return: The length of `vector`.
    :rtype: float
    """
    # Note: This seems faster than `numpy.linalg.norm`.
    return np.sqrt(np.dot(vector, vector))


def get_scale(left, right):
    """
    Calculate scaling factors based upon total counts. One of the factors
    is always one (the other is either one or larger than one).

    :arg array_like left, right: A vector.

    :return: A tuple of scaling factors.
    :rtype: float, float
    """
    left_scale = 1.0
    right_scale = 1.0

    # Calculate the scaling factors in such a way that no element is
    # between 0 and 1.
    left_sum = np.sum(left)
    right_sum = np.sum(right)

    if left_sum < right_sum:
        left_scale = right_sum / left_sum
    else:
        right_scale = left_sum / right_sum

    return left_scale, right_scale


def scale_down(left, right):
    """
    Normalise scaling factor between 0 and 1.

    :arg float left, right: Scaling factors.

    :return: Tuple of normalised scaling factors.
    :rtype: float, float
    """
    factor = max(left, right)

    return left / factor, right / factor


def positive(vector, mask):
    """
    Set all zero positions in `mask` to zero in `vector`.

    :arg array_like vector, mask: Vector.

    :return: `vector` with all zero positions in `mask` set to zero.
    :rtype: numpy.ndarray
    """
    return np.multiply(vector, np.asanyarray(mask, dtype=bool))


def multiset(left, right, pairwise):
    """
    Calculate the multiset distance between two vectors.

    :arg array_like left, right: Vector.
    :arg function pairwise: A pairwise distance function.

    :return: The multiset distance between `left` and `right`.
    :rtype: float

    Note that `function` must be vectorized, i.e., it is called directly on
    NumPy arrays, instead of on their pairwise elements. If your function only
    works on individual elements, convert it to a NumPy ufunc first. For
    example::

        >>> f = np.vectorize(f, otypes=['float'])
    """
    left = np.asanyarray(left)
    right = np.asanyarray(right)

    nonzero = np.where(np.logical_or(left, right))
    distances = pairwise(left[nonzero], right[nonzero])
    return distances.sum() / (len(distances) + 1)


def euclidean(left, right):
    """
    Calculate the Euclidean distance between two vectors.

    :arg array_like left, right: Vector.

    :return: The Euclidean distance between `left` and `right`.
    :rtype: float
    """
    return vector_length(np.subtract(left, right))


def cosine_similarity(left, right):
    """
    Calculate the Cosine similarity between two vectors.

    :arg array_like left, right: Vector.

    :return: The Cosine similarity between `left` and `right`.
    :rtype: float
    """
    return np.dot(left, right) / (vector_length(left) * vector_length(right))


#: Vector distance functions.
vector_distance = {
    "default": None,
    "euclidean": euclidean,
    "cosine": cosine_similarity
}


#: Pairwise distance functions. Arguments should be of type `numpy.ndarray`.
pairwise = {
    "prod": lambda x, y: abs(x - y) / ((x + 1) * (y + 1)),
    "sum": lambda x, y: abs(x - y) / (x + y + 1)
}


#: Summary functions.
summary = {
    "min": np.min,
    "average": np.mean,
    "median": np.median
}


#: Merge functions. Arguments should be of type `numpy.ndarray`.
mergers = {
    "sum": lambda x, y: x + y,
    "xor": lambda x, y: (x + y) * np.logical_xor(x, y),
    "int": lambda x, y: x * np.asanyarray(y, dtype=bool),
    "nint": lambda x, y: x * np.logical_not(y)
}
