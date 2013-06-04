#!/usr/bin/env python

#REQUIRES: MUSCLE

import os
import sys

lociFile = 'loci'

lociFile = open(lociFile, 'r')

for line in lociFile:

	locus = line.split('\t')
	print locus[1]

	variables = dict(
	fasta = locus[1] + '.fa', #your reference or reference sequences
	outfile = locus[1] + '.afa') #your outfile

	commands = """
	muscle3.8.31_i86win32.exe -in {fasta} -out {outfile}
	""".format(**variables)

	cmd_list = commands.split("\n")
	for cmd in cmd_list:
		os.system(cmd)
