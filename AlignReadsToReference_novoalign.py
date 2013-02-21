#!/usr/bin/env python

#REQUIRES: novoalign and samtools

#This script aligns your paired and unpaired reads to a reference using novoalign, and makes a pileup file using samtools

import os
import sys

#this is a wrap around for novoalign and samtools where each sample identifier was "index#" where # was a number between 1 - 50

ID = 'index' 
number = range(51)
number = number[1:51]

r1name = '_1p_final_renamed.fastq'  #extension of front reads
r2name = '_2p_final_renamed.fastq' #extension of back reads
uname = '_u_final_u_combined.fastq' #extension of unpaired reads

for num in number:

	variables = dict(
	reference = 'mitoRef.fa', #your reference or reference sequences
	#PARAMETERS for NOVOALIGN
	insertSize='160', #mean fragment length of your libraries
	insertSTDEV = '60', #standard deviation of length of your libraries
	aScore = '690', #stringency score...I don't remember what 90 means, but it is stringent
	a1 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC', #adapter sequence
	a2 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA', #adapter sequence
	truncatePaired ='100',
	truncateUnpaired ='200',
	read1 = ID+ str(num) + r1name,
	read2 = ID+ str(num) + r2name,
	unpaired = ID+str(num)+uname,
	novoOutPaired = ID+str(num) + '_out_paired',
	novoOutUnpaired = ID+str(num) + '_out_unpaired',
	outfile = ID+str(num)+'_mito_sorted') #name your output

#in order,
#1. index your reference sequence
#2. align paired reads to reference
#3. align unpaired reads to reference
#4. remove paired reads not aligned to reference
#5. remove unpaired reads not aligned to reference
#6. make a bam file for paired reads
#7. make a bam file for unpaired reads
#8. merge paired and unpaired reads
#9. sort the sam file
#10. 'index the file'
#11. make a pileup file.

	commands = """
	novoindex {reference}.ndx {reference}
	novoalign -d {reference}.ndx -f {read1} {read2} -i PE {insertSize} {insertSTDEV} -t {aScore} -n {truncatePaired} -a {a1} {a2} -F STDFQ -o SAM > {novoOutPaired}
	novoalign -d {reference}.ndx -f {unpaired} -t {aScore} -n {truncateUnpaired} -a {a1} {a2} -F STDFQ -o SAM > {novoOutUnpaired}
	grep -v ZS:Z:NM {novoOutPaired} > {novoOutPaired}.sam
	grep -v ZS:Z:NM {novoOutUnpaired} > {novoOutUnpaired}.sam
	samtools view -bS {novoOutPaired}.sam > {novoOutPaired}.bam
	samtools view -bS {novoOutUnpaired}.sam > {novoOutUnpaired}.bam
	samtools merge -f raw.bam {novoOutPaired}.bam {novoOutUnpaired}.bam
	samtools sort raw.bam {outfile}
	samtools index {outfile}.bam
	samtools mpileup -f {reference} {outfile}.bam > {outfile}.pileup
	""".format(**variables)

	cmd_list = commands.split("\n")
	for cmd in cmd_list:
		os.system(cmd)

	os.system('rm raw.bam *.ndx') #*out_paired* *out_unpaired*










