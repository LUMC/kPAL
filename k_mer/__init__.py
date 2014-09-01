"""
k_mer: Analysis toolkit and programming library for k-mer profiles.


Copyright (c) 2013 Leiden University Medical Center <humgen@lumc.nl>
Copyright (c) 2013 Jeroen F.J. Laros <j.f.j.laros@lumc.nl>

Licensed under the MIT license, see the LICENSE file.
"""

import argparse
import os

import h5py

__version_info__ = ('0', '3', '0')

__version__ = '.'.join(__version_info__)
__author__ = 'LUMC, Jeroen F.J. Laros'
__contact__ = 'J.F.J.Laros@lumc.nl'
__homepage__ = 'https://git.lumc.nl/j.f.j.laros/k-mer'

usage = __doc__.split("\n\n\n")

class ProtectedFileType(argparse.FileType):
    def __call__(self, string):
        if 'w' in self._mode and os.path.exists(string):
            raise IOError('failed to create "{}": file exists.'.format(string))
        return super(ProtectedFileType, self).__call__(string)

class ProfileFileType(object):
    """Factory for creating file object types

    Instances of FileType are typically passed as type= arguments to the
    ArgumentParser add_argument() method.

    Keyword Arguments:
        - mode -- A string indicating how the file is to be opened. Accepts the
            same values as the builtin open() function.

    Todo: Tidy and docstring.
    """
    def __init__(self, mode='r'):
        self._mode = mode

    def __call__(self, string):
        if 'w' in self._mode and os.path.exists(string):
            raise IOError('failed to create "{}": file exists.'.format(string))
        try:
            handle = h5py.File(string, self._mode)
            if 'w' in self._mode:
                handle.attrs['format'] = 'kMer'
                handle.attrs['version'] = '1.0.0'  # Todo: module variable
                handle.attrs['producer'] = 'kMer %s' % __version__
                handle.create_group('profiles')
            # Todo: When reading and/or appending, check if format is 'kMer'
            # and version is supported.
            return handle
        except IOError as e:
            message = "can't open '%s': %s"
            raise ArgumentTypeError(message % (string, e))

    def __repr__(self):
        args = self._mode, self._bufsize
        args_str = ', '.join(repr(arg) for arg in args if arg != -1)
        return '%s(%s)' % (type(self).__name__, args_str)

def doc_split(func):
    return func.__doc__.split("\n\n")[0]

def version(name):
    return "%s version %s\n\nAuthor   : %s <%s>\nHomepage : %s" % (name,
        __version__, __author__, __contact__, __homepage__)
