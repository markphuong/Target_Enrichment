#!/usr/bin/env python

import dendropy
import os
import sys

ID = 'index' 
number = range(45)
number = number[1:45]

lociFile = 'loci'

lociFile = open(lociFile, 'r')

for line in lociFile:

	locus = line.split('\t')
	locusOutFile = locus[1] + '_noHistoric.fa'
	locusOutFile = open(locusOutFile, 'w')

	locusInFile = locus[1] +'.fa'

	seqFile = dendropy.DnaCharacterMatrix.get_from_path(locusInFile, "fasta")

	for num in number:
		seqFile
