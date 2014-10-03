"""
kMer: Analysis toolkit and programming library for k-mer profiles.


Copyright (c) 2013-2014 Leiden University Medical Center <humgen@lumc.nl>
Copyright (c) 2013-2014 Jeroen F.J. Laros <j.f.j.laros@lumc.nl>
Copyright (c) 2014 Martijn Vermaat <m.vermaat.hg@lumc.nl>

Licensed under the MIT license, see the LICENSE file.
"""


from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import str
from future.utils import native

import argparse
from io import open
import os
import sys

import h5py
import semantic_version


__version_info__ = ('1', '0', '1')
__date__ = '3 Oct 2014'


__version__ = '.'.join(__version_info__)
__author__ = 'LUMC, Jeroen F.J. Laros, Martijn Vermaat'
__contact__ = 'J.F.J.Laros@lumc.nl'
__homepage__ = 'https://git.lumc.nl/j.f.j.laros/k-mer'


USAGE = __doc__.split("\n\n\n")
FORMAT_VERSION = semantic_version.Version('1.0.0')
FORMAT_ACCEPT = semantic_version.Spec('>=1.0.0,<2.0.0')


# Same as argparse.FileType in Python 3, but using io.open to get the same
# behaviour on Python 2, and adding protection against overwriting existing
# files.
class FileType(object):
    def __init__(self, mode='r', bufsize=-1, encoding=None, errors=None):
        self._mode = mode
        self._bufsize = bufsize
        self._encoding = encoding
        self._errors = errors

    def __call__(self, string):
        # the special argument "-" means sys.std{in,out}
        if string == '-':
            if 'r' in self._mode:
                return sys.stdin
            elif 'w' in self._mode:
                return sys.stdout
            else:
                msg = 'argument "-" with mode %r' % self._mode
                raise ValueError(msg)

        # all other arguments are used as file names
        try:
            if 'w' in self._mode and os.path.exists(string):
                raise OSError('file exists')
            return open(string, self._mode, self._bufsize, self._encoding,
                        self._errors)
        except OSError as e:
            message = "can't open '%s': %s"
            raise argparse.ArgumentTypeError(message % (string, e))

    def __repr__(self):
        args = self._mode, self._bufsize
        kwargs = [('encoding', self._encoding), ('errors', self._errors)]
        args_str = ', '.join([repr(arg) for arg in args if arg != -1] +
                             ['%s=%r' % (kw, arg) for kw, arg in kwargs
                              if arg is not None])
        return '%s(%s)' % (type(self).__name__, args_str)


class ProfileFileType(object):
    def __init__(self, mode='r'):
        self._mode = mode

    def __call__(self, string):
        try:
            if 'w' in self._mode and os.path.exists(string):
                raise IOError('file exists')
            handle = h5py.File(string, self._mode)
            if 'w' in self._mode:
                handle.attrs['format'] = 'kMer'
                handle.attrs['version'] = native(str(FORMAT_VERSION))
                handle.attrs['producer'] = 'kMer %s' % __version__
                handle.create_group('profiles')
            else:
                if handle.attrs.get('format') != 'kMer':
                    raise IOError('not a k-mer profile file')
                version = semantic_version.Version(handle.attrs['version'])
                if version not in FORMAT_ACCEPT:
                    raise IOError('file format version %s not supported'
                                  % version)
            return handle
        except IOError as e:
            message = "can't open '%s': %s"
            raise argparse.ArgumentTypeError(message % (string, e))

    def __repr__(self):
        return '%s(%s)' % (type(self).__name__, self._mode)


def doc_split(func):
    return func.__doc__.split("\n\n")[0]


def version(name):
    return '{0} version {1}\n\nAuthor   : {2} <{3}>\nHomepage : {4}'.format(
        name, __version__, __author__, __contact__, __homepage__)
