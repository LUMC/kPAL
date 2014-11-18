from setuptools import setup

try:
    with open('README.rst') as readme:
        long_description = readme.read()
except IOError:
    long_description = 'kMer has been renamed to *kPAL*'

setup(
    name='kMer',
    version='1.0.2',
    description='k-mer analysis toolkit and programming library.',
    long_description=long_description,
    author='LUMC',
    author_email='j.f.j.laros@lumc.nl',
    url='https://github.com/LUMC/kPAL',
    license='MIT License',
    platforms=['any'],
    packages=['k_mer'],
    install_requires=['kPAL'],
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
