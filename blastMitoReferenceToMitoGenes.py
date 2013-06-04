#!/usr/bin/env python

#REQUIRES: blastn

import os
import sys

ID = 'index' 
number = range(51)
number = number[1:51]

for num in number:

	sample = ID + str(num)
	print sample + 'LESGO'

	variables = dict(
	reference = 'mitoLoci100plus.fasta', #your reference or reference sequences
	infile = sample + '_finalMitoContig.fa', #your file extensions from your identifier (here, the identifier is Index, and the extension is _mitoseq.fa
	evalue = 10000, #stringency threshold for blastn
	wordsize = 22, #helps optimize blast finding short sequences
	outfile = ID + str(num) + '_mitoBlastOutput') #your outfile

	commands = """
	blastn -db {reference} -query {infile} -outfmt 6 -evalue {evalue} -out {outfile} -word_size {wordsize}
	""".format(**variables)

	cmd_list = commands.split("\n")
	for cmd in cmd_list:
		os.system(cmd)
