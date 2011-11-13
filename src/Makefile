# Makefile
#
# qmake -v PATH -cwd -- -j 100
#

# SGE configuration.
SGE_PE = BWA

# k-mer tools.
KMER_DIR = /home/jfjlaros/projects/k-mer/trunk/src
KMER = python $(KMER_DIR)/kMer.py
KMERMERGE = python $(KMER_DIR)/kMerMerge.py
KMERMATRIX = python $(KMER_DIR)/kMerMatrix.py
KMER_SIZE = 9
KMER_PRECISION = 10

# Find the input files.
INPUT = $(shell ls *.fq)
SAMPLES = $(shell ls *_1.fq | sed 's/_1.fq//')


all: matrix.$(KMER_SIZE)

%.count: %.fq
	$(KMER) -i $< -o $@ -k $(KMER_SIZE)

%.merged: %_1.count %_2.count
	$(KMERMERGE) -i $^ -o $@

matrix.$(KMER_SIZE): $(addsuffix .merged, $(SAMPLES))
	$(KMERMATRIX) -i $^ -o $@ -p $(KMER_PRECISION)
