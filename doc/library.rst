Using the Python library
========================

kMer provides a light-weight Python library for creating, analysing, and
manipulating *k*-mer profiles. It is implemented on top of `NumPy
<http://www.numpy.org/>`_.

This is a gentle introduction to the library. Consult the :ref:`api` for more
detailed documentation.


*k*-mer profiles
----------------

.. currentmodule:: k_mer.klib

The class :class:`Profile` is the central object in kMer. It encapsulates
*k*-mer counts and provides operations on them.

Instead of using the :class:`Profile` constructor directly, you should
generally use one of the profile construction methods. One of those is
:meth:`Profile.from_fasta`. The following code creates a 6-mer profile by
counting from a FASTA file::

    >>> from k_mer.klib import Profile
    >>> p = Profile.from_fasta(open('a.fasta'), 6)

The profile object has several properties. For example, we can ask for the
*k*-mer length (also known as *k*), the total *k*-mer count, or the median
count per *k*-mer::

    >>> p.length
    6
    >>> p.total
    49995
    >>> p.median
    12.0

Counts are stored as a NumPy :class:`~numpy.ndarray` of integers, one for each
possible *k*-mer, in alphabetical order::

    >>> len(p.counts)
    4096
    >>> p.counts
    array([ 8, 11,  5, ...,  7, 12, 13])

We can get the index in that array for a certain *k*-mer using the
:meth:`~Profile.dna_to_binary` method::

    >>> i = p.dna_to_binary('AATTAA')
    >>> p.counts[i]
    13


Storing *k*-mer profiles
------------------------

Todo.

.. Todo: Explain ProfileFileType.

    We can save the profile to a file (see :ref:`fileformat`) for later use::


        >>> from k_mer import ProfileFileType
        >>> with ProfileFileType('w')('a.k6') as f:
        ...     p.save(f)

    And read it::

        >>> with ProfileFileType()('a.k6') as f:
        ...     p = Profile.from_file(f)


Differences between *k*-mer profiles
------------------------------------

Todo.
