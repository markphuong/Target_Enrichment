#!/usr/bin/env python

#REQUIRES: raxmlHPC


import os
import sysz

lociFile = 'loci'

lociFile = open(lociFile, 'r')

for line in lociFile:

	locus = line.split('\t')

	variables = dict(
	phylip = locus[1] + '.phy', #the phylip file you want to put in
	marker = locus[1], #the locus you are on
	outgroup = 'index31_' + locus[1] + '_variegatus,' + 'index9_' + locus[1] +'_variegatus,' + 'index32_' + locus[1] + '_variegatus') #your outgroup

	commands = """
	raxmlHPC -f a -s {phylip} -x 12345 -# 100 -m GTRCAT -n {marker} -o {outgroup}
	""".format(**variables)

	cmd_list = commands.split("\n")
	for cmd in cmd_list:
		os.system(cmd)

