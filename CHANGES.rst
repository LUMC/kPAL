Changelog
=========

This is a record of changes made between each kMer release.


Version 1.0.0
-------------

Released on October 2nd, 2014.

- Also count *k*-mers if *k* equals the length of the string.
- Python 2.6 compatibility.
- Added unit tests and a `tox <https://testrun.org/tox/>`_ configuration.
- Use a NumPy ndarray for storing *k*-mer counts.
- New multi-profile HDF5 file format (see :ref:`fileformat`).
- Fix splitting a profile for calculating balance. Palindromes were previously
  not taken into account when splitting a profile. We now double all counts,
  so palindrome counts can be evenly distributed over both sides (see `GitLab
  #1 <https://git.lumc.nl/j.f.j.laros/k-mer/issues/1>`_).
- Our own implementation of a vector's median contained two bugs. Better to
  use a library for this.
- Fix Euclidean distance between two vectors. Don't add one to the sum of
  squares. The distance between two empty vectors should be 0, not 1.
- Rename `k_mer.klib.kMer` to `k_mer.klib.Profile`.
- Support Python 3.3 and 3.4 (:ref:`see below <v1.0.0-py3>`).
- Generalize custom function arguments (:ref:`see below <v1.0.0-custom>`).
- `Travis CI <https://travis-ci.org/LUMC/kMer>`_ configuration.
- `Sphinx documentation <http://kmer.readthedocs.org/>`_ including a user
  guide and API reference.
- Renamed the `index` command to `count`.
- Renamed the `diff` command to `distance` and the builtin pairwise distance
  functions from `diff-prod` and `diff-sum` to just `prod` and `sum`.


.. _v1.0.0-py3:

Support Python 3.3 and 3.4
^^^^^^^^^^^^^^^^^^^^^^^^^^

*TL;DR:* kMer supports Python 2 and 3 and every module has the following line
at the top::

    >>> from __future__ import (absolute_import, division, print_function,
                                unicode_literals)

We now support Python versions 2.6, 2.7, 3.3, and 3.4 in a single codebase
without using 2to3. We don't support Python 3.2 because BioPython does not.

We use the `Python future <http://python-future.org/>`_ package as a
compatibility layer between Python 2 and Python 3. The goal is to use a
single, clean Python 3.x-compatible codebase to support both Python 2 and
Python 3 with minimal overhead.

Most changes are quite straightforward (e.g., absolute imports, print
statement, division operator). The main painpoint is of course the bytestring
versus unicode story. We now `import unicode_literals
<http://python-future.org/imports.html#should-i-import-unicode-literals>`_ in
each module and maintain that all text in kMer is unicode (`unicode` in Python
2, `str` in Python 3).


.. _v1.0.0-custom:

Generalize custom function arguments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Custom function arguments in the command line interface can now be either a
Python expression or importable name. For example, all commands accepting a
summary function argument, also except a custom summary function argument
which should be one of:

1. A Python expression over the NumPy ndarray `values` (e.g.,
   ``np.max(values)``).
2. An importable name (e.g., ``package.module.summary``) that can be called
   with an ndarray as argument.

Likewise for custom merger and pairwise functions (here the expression is over
the two NumPy ndarrays `left` and `right`).


Version 0.3.0
-------------

Released on July 3rd, 2014.

- Usage of the Euclidean distance is now handled differently, breaking
  backwards compatibility.
- Added Cosine similarity measure and generalised distance parameters.
- Fixed broken setup script.
- Added custom merging functionality.


Version 0.2.0
-------------

Released on March 23rd, 2014.

- New command line interface, using positional arguments for required
  parameters.
- Added checking for existing files to prevent overwriting them.
- Fixed a bug in the scale subcommand that prevented scaling.
- Added a version parameter.
- Updated the homepage.
- Made code PEP 8 compliant.
- Switched to Sphynx docstrings.
- Added keyword selection for distance and smoothing functions.
- Added support for custom distance and smoothing functions.
- Added CHANGELOG and README.


Version 0.1.0
-------------

Released on September 24th, 2013.

- Start of log.
