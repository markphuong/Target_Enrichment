#!/usr/bin/env python

#REQUIRES: blastn

import os
import sys
import dendropy

lociFile = dendropy.DnaCharacterMatrix.get_from_path("mitoLoci.fa", "fasta")

for loci in lociFile:

	loci = str(loci)
#	sample = ID + str(num)
	print loci + ' LESGO'

	variables = dict(
	reference = 'mitoLoci.fa', #your reference or reference sequences
	infile = loci + '.fa', #your file extensions from your identifier (here, the identifier is Index, and the extension is _mitoseq.fa
	evalue = 10000, #stringency threshold for blastn
	wordsize = 22, #helps optimize blast finding short sequences
	outfile = loci + '_mitoBlastOutput') #your outfile

	commands = """
	blastn -db {reference} -query {infile} -outfmt 6 -evalue {evalue} -out {outfile} -word_size {wordsize}
	""".format(**variables)

	cmd_list = commands.split("\n")
	for cmd in cmd_list:
		os.system(cmd)
