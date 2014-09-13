"""
k_mer: Analysis toolkit and programming library for k-mer profiles.


Copyright (c) 2013 Leiden University Medical Center <humgen@lumc.nl>
Copyright (c) 2013 Jeroen F.J. Laros <j.f.j.laros@lumc.nl>

Licensed under the MIT license, see the LICENSE file.
"""

import argparse
import os

import h5py
import semantic_version

__version_info__ = ('0', '3', '0')

__version__ = '.'.join(__version_info__)
__author__ = 'LUMC, Jeroen F.J. Laros'
__contact__ = 'J.F.J.Laros@lumc.nl'
__homepage__ = 'https://git.lumc.nl/j.f.j.laros/k-mer'

USAGE = __doc__.split("\n\n\n")
FORMAT_VERSION = semantic_version.Version('1.0.0')
FORMAT_ACCEPT = semantic_version.Spec('>=1.0.0,<2.0.0')

class ProtectedFileType(argparse.FileType):
    def __call__(self, string):
        if 'w' in self._mode and os.path.exists(string):
            raise IOError('failed to create "{}": file exists.'.format(string))
        return super(ProtectedFileType, self).__call__(string)

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
                handle.attrs['version'] = str(FORMAT_VERSION)
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
    return "%s version %s\n\nAuthor   : %s <%s>\nHomepage : %s" % (name,
        __version__, __author__, __contact__, __homepage__)
