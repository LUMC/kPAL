import os
from setuptools import setup
import sys

if sys.version_info < (2, 6):
    raise Exception('kMer requires Python 2.6 or higher.')

dependencies = ['numpy', 'h5py', 'biopython', 'semantic-version', 'future']
if sys.version_info[:2] == (2, 6):
    dependencies.extend(['argparse', 'importlib'])

try:
    with open('README.md') as readme:
        long_description = readme.read()
except IOError:
    long_description = 'See https://pypi.python.org/pypi/kmer'

# This is quite the hack, but we don't want to import our package from here
# since that's recipe for disaster (it might have some uninstalled
# dependencies, or we might import another already installed version).
distmeta = {}
for line in open(os.path.join('k_mer', '__init__.py')):
    try:
        field, value = (x.strip() for x in line.split('='))
    except ValueError:
        continue
    if field == '__version_info__':
        value = value.strip('[]()')
        value = '.'.join(x.strip(' \'"') for x in value.split(','))
    else:
        value = value.strip('\'"')
    distmeta[field] = value

setup(
    name='kMer',
    version=distmeta['__version_info__'],
    description='k-mer analysis toolkit and programming library.',
    long_description=long_description,
    author=distmeta['__author__'],
    author_email=distmeta['__contact__'],
    url=distmeta['__homepage__'],
    license='MIT License',
    platforms=['any'],
    packages=['k_mer'],
    install_requires=dependencies,
    entry_points={
        'console_scripts': ['kMer = k_mer.kmer:main']
        },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    keywords='bioinformatics'
)
