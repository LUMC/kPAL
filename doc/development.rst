.. highlight:: none

Development
===========

Development of kPAL happens on GitHub:
https://github.com/LUMC/kPAL


Contributing
------------

Contributions to kPAL are very welcome! They can be feature requests, bug
reports, bug fixes, unit tests, documentation updates, or anything els you may
come up with.

Start by installing all kPAL development dependencies::

    $ pip install -r requirements.txt

This installs dependencies for building the documentation and running unit
tests.

After that you'll want to install kPAL in *development mode*::

    $ pip install -e .

.. note:: Instead of copying the source code to the installation directory,
          this only links from the installation directory to the source code
          such that any changes you make to it are directly available in the
          environment.


Documentation
-------------

The `latest documentation <http://kpal.readthedocs.org/>`_ with user guide and
API reference is hosted at Read The Docs.

You can also compile the documentation directly from the source code by
running ``make html`` from the ``doc/`` subdirectory. This requires `Sphinx`_
to be installed.


Unit tests
----------

To run the unit tests with `pytest`_, just run::

    $ py.test

Use `tox`_ to run the unit tests in all supported Python environments
automatically::

    $ tox


Coding style
------------

In general, try to follow the `PEP 8`_ guidelines for Python code and `PEP
257`_ for docstrings.

You can use the `flake8`_ tool to assist in style and error checking.


Versioning
----------

A normal version number takes the form X.Y.Z where X is the major version, Y
is the minor version, and Z is the patch version. Development versions take
the form X.Y.Z.dev where X.Y.Z is the closest future release version.

Note that this scheme is not 100% compatible with `SemVer`_ which would
require X.Y.Z-dev instead of X.Y.Z.dev but `compatibility with setuptools
<http://peak.telecommunity.com/DevCenter/setuptools#specifying-your-project-s-version>`_
is more important for us. Other than that, version semantics are as described
by SemVer.

Releases are `published at PyPI <https://pypi.python.org/pypi/kPAL>`_ and
available from the git repository as tags.


Release procedure
^^^^^^^^^^^^^^^^^

Releasing a new version is done as follows:

1. Make sure the section in the ``CHANGES.rst`` file for this release is
   complete and there are no uncommitted changes.

   .. note::

    Commits since release X.Y.Z can be listed with ``git log vX.Y.Z..`` for
    quick inspection.

2. Update the ``CHANGES.rst`` file to state the current date for this release
   and edit ``kpal/__init__.py`` by updating `__date__` and removing the
   ``dev`` value from `__version_info__`.

   Commit and tag the version update::

       git commit -am 'Bump version to X.Y.Z'
       git tag -a 'vX.Y.Z'
       git push --tags

3. Upload the package to PyPI::

       python setup.py sdist upload

4. Add a new entry at the top of the ``CHANGES.rst`` file like this::

       Version X.Y.Z+1
       ---------------

       Release date to be decided.

   Increment the patch version and add a ``dev`` value to `__version_info__`
   in ``kpal/__init__.py`` and commit these changes::

       git commit -am 'Open development for X.Y.Z+1'


.. _Sphinx: http://sphinx-doc.org/
.. _pytest: http://pytest.org/
.. _tox: https://testrun.org/
.. _PEP 8: http://www.python.org/dev/peps/pep-0008/
.. _PEP 257: http://www.python.org/dev/peps/pep-0257/
.. _flake8: http://flake8.readthedocs.org/en/latest/
.. _SemVer: http://semver.org/
