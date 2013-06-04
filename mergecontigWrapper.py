#!/usr/bin/env python

import os
import sys

ID = 'index' 
number = range(51)
number = number[1:51]

for num in number:

	variables = dict(
	sample = ID + str(num)+'_contig.fa')

	commands = """
	perl 6finalAssemblyCommented.pl -a {sample}
	""".format(**variables)

	cmd_list = commands.split("\n")
	for cmd in cmd_list:
		os.system(cmd)
