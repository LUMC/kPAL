.. _fileformat:

*k*-mer profile file format
===========================

The file format kPAL uses to store *k*-mer profiles is `HDF5
<http://www.hdfgroup.org/>`_. Here we describe the structure within a *k*-mer
profile file.


Versioning
----------

The file format is versioned roughly according to `semantic versioning
<http://semver.org/>`_. Software designed to work with files in version
*MAJOR.MINOR.PATCH* should be able to work with files in later versions with
the same *MAJOR* version without modification.


Current version: 1.0.0
----------------------

The HDF5 toplevel attributes are:

- **format** (`string`) -- This is always set to ``kMer``.
- **version** (`string`) -- Currently ``1.0.0``.
- **producer** (`string`) -- Anything, for example ``My k-mer program 1.2.1``.

Each *k*-mer profile is a dataset under the ``/profiles`` group, named
``/profiles/<profile_name>``. The data is a one-dimensional array of integers
of length :math:`4^k` (where :math:`k` is the *k*-mer length) and is gzip
compressed. This dataset has the following attributes:

- **length** (`integer`): *k*-mer length (also know as *k*).
- **total** (`integer`): Sum of *k*-mer counts.
- **non_zero** (`integer`): Number of *k*-mers with a non-zero count.
- **mean** (`float`): Mean of *k*-mer counts.
- **median** (`integer`): Median of *k*-mer counts.
- **std** (`float`): Standard deviation of *k*-mer counts.

Within one file, all profiles must have the same value for the `length`
attribute.

All strings and object names in the file are unicode strings encoded as
described in the `h5py documentation
<http://docs.h5py.org/en/latest/strings.html>`_.


Changes from older versions
---------------------------

None yet.
