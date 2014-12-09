.. highlight:: none

.. _install:

Installation
============

The kPAL source code is `hosted on GitHub
<https://github.com/LUMC/kPAL>`_. Supported Python versions for running kPAL
are 2.6, 2.7, 3.3, and 3.4. kPAL can be installed either via the Python
Package Index (PyPI) or from the source code.


Dependencies
------------

kPAL depends on the following Python libraries:

- `NumPy <http://www.numpy.org/>`_
- `h5py <http://www.h5py.org/>`_
- `biopython <http://biopython.org/>`_

The easiest way to use kPAL is with the `Anaconda distribution
<https://store.continuum.io/cshop/anaconda/>`_ which comes with these
libraries installed.

Alternatively, you can install them using their binary packages for your
operating system.

Although all dependencies will also be automatically installed if they aren't
yet when installing kPAL, you may still want to have them installed
beforehand. Automatic installation requires compilation from source, which
takes a lot of time and needs several compilers and development libraries to
be available. The options noted above are often much more convenient.


Latest kPAL release
-------------------

To install the latest release from `PyPI <https://pypi.python.org/pypi/kPAL>`_
using pip::

    pip install kPAL


kPAL development version
------------------------

You can also clone and use the latest development version directly from the
GitHub repository::

    git clone https://github.com/LUMC/kPAL.git
    cd kPAL
    pip install -e .
