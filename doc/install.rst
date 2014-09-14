.. highlight:: none

.. _install:

Installation
============

The kMer source code is `hosted on GitLab
<https://git.lumc.nl/j.f.j.laros/k-mer>`_. Supported Python versions for
running kMer are 2.6, 2.7, 3.3, and 3.4. kMer can be installed either via the
Python Package Index (PyPI) or from the source code.


Dependencies
------------

kMer depends on the following Python libraries:

- `NumPy <http://www.numpy.org/>`_
- `h5py <http://www.h5py.org/>`_
- `biopython <http://biopython.org/>`_

The easiest way to use kMer is with the `Anaconda distribution
<https://store.continuum.io/cshop/anaconda/>`_ which comes with these
libraries installed.

Alternatively, you can install them using their binary packages for your
operating system.

Although all dependencies will also be automatically installed if they aren't
yet when installing kMer, you may still want to have them installed
beforehand. Automatic installation requires compilation from source, which
takes a lot of time and needs several compilers and development libraries to
be available. The options noted above are often much more convenient.


Latest kMer release
-------------------

To install the latest release from `PyPI <https://pypi.python.org/pypi/kMer>`_
using pip::

    pip install kMer


kMer development version
------------------------

You can also clone and use the latest development version directly from the
GitLab repository::

    git clone https://git.lumc.nl/j.f.j.laros/k-mer.git
    cd k-mer
    pip install -e .
