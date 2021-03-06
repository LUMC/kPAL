.. highlight:: none

.. _tutorial:

Tutorial
========

Before following this tutorial, make sure kPAL is installed  properly::

    $ kpal -h

This should print a help message. If it does not, follow :ref:`install`.

We work with an artificial dataset consisting of 200 *read pairs* from four
different samples. They are randomly generated so have no biological
relevance.

.. note:: Download the data: :download:`tutorial.zip <downloads/tutorial.zip>`

    Now unzip the file and go to the resulting directory::

        $ unzip -q tutorial.zip
        $ cd tutorial
        $ ls
        a_1.fa  a_2.fa  b_1.fa  b_2.fa  c_1.fa  c_2.fa  d_1.fa  d_2.fa

We'll create *k*-mer profiles for these samples and try to compare them.


*k*-mer counting
----------------

kPAL can count *k*-mers in any number of fasta files and store the results in
one *k*-mer profile file. By default, the profiles in the file are named
according to the original fasta filenames.

Let's count 8-mers in the first read for all samples and write the profiles to
``reads_1.k8``::

    $ kpal count -k 8 *_1.fa reads_1.k8

Using the `info` command, we can get an overview of our profiles::

    $ kpal info reads_1.k8
    File format version: 1.0.0
    Produced by: kPAL 2.0.0

    Profile: a_1
    - k-mer length: 8 (65536 k-mers)
    - Zero counts: 49395
    - Non-zero counts: 16141
    - Sum of counts: 18600
    - Mean of counts: 0.284
    - Median of counts: 0.000
    - Standard deviation of counts: 0.535

    Profile: b_1
    - k-mer length: 8 (65536 k-mers)
    - Zero counts: 49348
    - Non-zero counts: 16188
    - Sum of counts: 18600
    - Mean of counts: 0.284
    - Median of counts: 0.000
    - Standard deviation of counts: 0.533

    Profile: c_1
    - k-mer length: 8 (65536 k-mers)
    - Zero counts: 49388
    - Non-zero counts: 16148
    - Sum of counts: 18600
    - Mean of counts: 0.284
    - Median of counts: 0.000
    - Standard deviation of counts: 0.534

    Profile: d_1
    - k-mer length: 8 (65536 k-mers)
    - Zero counts: 49345
    - Non-zero counts: 16191
    - Sum of counts: 18600
    - Mean of counts: 0.284
    - Median of counts: 0.000
    - Standard deviation of counts: 0.533


Merging profiles
----------------

For completeness, we also want to include *k*-mer counts for the second read
in our analysis. We can do so using the `merge` command::

    $ kpal count -k 8 *_2.fa reads_2.k8
    $ kpal merge reads_1.k8 reads_2.k8 merged.k8

.. note:: Merging two *k*-mer profiles this way is equivalent to first
          concatenating both fasta files and counting in the result.

By default, profiles from both files are merged pairwise in alphabetical
order. If you need another pairing, you can provide profile names to use for
both files. For example, the following is a more explicit version of the
previous command::

    $ kpal merge reads_1.k8 reads_2.k8 merged.k8 -l a_1 b_1 c_1 d_1 -r a_2 b_2 c_2 d_2

We can check that, indeed, the total *k*-mer count has doubled compared to our
previous numbers::

    $ kpal info merged.k8 -p c_1_c_2
    File format version: 1.0.0
    Produced by: kPAL 2.0.0

    Profile: c_1_c_2
    - k-mer length: 8 (65536 k-mers)
    - Zero counts: 37138
    - Non-zero counts: 28398
    - Sum of counts: 37200
    - Mean of counts: 0.568
    - Median of counts: 0.000
    - Standard deviation of counts: 0.753


Distance between profiles
-------------------------

We can compare two profiles by using a distance function. By default,
`distance` uses the multiset distance parameterised by the `prod` pairwise
distance function (:math:`f_2` in :ref:`method-distance`)::

    $ kpal distance reads_1.k8 reads_2.k8 -l c_1 -r c_2
    c_1 c_2 0.456

All profiles in a file can be compared pairwise to produce a distance matrix
with the `matrix` command. It first writes the number of profiles compared
followed by their names, and then the distance matrix itself. Here we ask it
to print the result to standard output (using ``-`` for the output filename)::

    $ kpal matrix merged.k8 -
    4
    a_1_a_2
    b_1_b_2
    c_1_c_2
    d_1_d_2
    0.415
    0.416 0.416
    0.414 0.413 0.414


Enforcing strand balance
------------------------

Todo.


Custom merge functions
----------------------

Todo.
