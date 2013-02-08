#!/usr/bin/env python



import sys
import os
import subprocess as sp
import shlex

ID = 'index'
number = range(50)
number = number[1:51]

for num in number:
	
	sample = ID + str(num)
	outpaired =  sample+'_out_paired.sam'
	outunpaired = sample + '_out_unpaired.sam'
	read1 = sample+'_1p_final_renamed.fastq'
	read2 = sample+'_2p_final_renamed.fastq'
	unpaired = sample + '_u_final_u_combined.fastq'

	variables = [outpaired, outunpaired, read1, read2, unpaired]

	values = []
	for var in variables:
		cmd = "grep HS1 "+ var+" | wc -l"
		p = sp.Popen(cmd, stdout=sp.PIPE, shell=True)
		output = p.communicate()[0]
		values.append(output.strip("\n"))

	values = map(float, values)

	specificity = (values[0] + values[1]) / ( values[2] + values[3] + values[4])

	out = open('mito_specificity.txt', 'a')
	out.write(ID, '\t', str(specificity))
	out.close()
