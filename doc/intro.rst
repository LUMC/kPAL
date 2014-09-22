.. highlight:: none

Introduction
============

kMer provides a command-line analysis toolkit for creating, analysing, and
manipulating *k*-mer profiles. It is implemented in Python.

After following :ref:`install`, kMer can be started by typing::

    $ kMer

More information about the available commands and their arguments is printed
by adding the `-h` argument.

For example, to count all 9-mers in a FASTA file, use the `count` command::

    $ kMer count -k 9 example.k9 example.fasta

Below, we provide an overview of all functions of kMer that are available via
the command-line interface:

===========  =================================================================
Command      Description
===========  =================================================================
count        Make a profile from a FASTA file.
merge        Merge two profiles.
balance      Balance a profile on the frequency of *k*-mers and their reverse
             complements.
showbalance  Calculate the balance of a profile.
meanstd      Show the mean and standard deviation of *k*-mer frequencies.
distr        Calculate the distribution of the frequencies in a profile.
info         Print basic statistics on a given profile.
getcount     Retrieve the count for a particular *k*-mer.
positive     Only keep counts that are positive in both profiles.
scale        Scale profiles such that the total number of *k*-mer frequencies
             is equal.
shrink       Shrink a profile, effectively reducing *k*-mer length.
shuffle      Randomise a profile.
smooth       Smooth two profiles by collapsing sub-profiles.
distance     Calculate the distance between two profiles.
matrix       Make a pairwise distance matrix for a series of *k*-mer profiles.
===========  =================================================================

More information about the methods implemented by kMer can be found in
:ref:`method`. Some examples of working with the toolkit are shown in
:ref:`tutorial`.
