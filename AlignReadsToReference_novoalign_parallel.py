#!/usr/bin/env python

#REQUIRES: novoalign and samtools
#REQUIRES: a map file, with first column as sample ID, and second file as which fasta it goes to. The reason you have different fastas for different samples is because of divergent mtDNA genomes
#elements in the map file are separated by a tab

#This script aligns your paired and unpaired reads to a reference using novoalign, and makes a pileup file using samtools

import os
import sys
import argparse
import multiprocessing

#this is a wrap around for novoalign and samtools where each sample identifier was "index#" where # was a number between 1 - 50

def get_args(): #arguments needed to give to this script
	parser = argparse.ArgumentParser(description="run novoalign")

	#forces required argument to let it run
	required = parser.add_argument_group("required arguments") 
	required.add_argument("--map", help="textfile with samples to run and what fasta file to match it to", required=True) #A map file with the sample ID and the fasta file it goes to

	return parser.parse_args()

def align(element):

	ID = element[0]
	fastafile = element[1]
	
	r1name = '_1p_final_renamed.fastq'  #extension of front reads
	r2name = '_2p_final_renamed.fastq' #extension of back reads
	uname = '_u_final_u_combined.fastq' #extension of unpaired reads

	variables = dict(
	sample = ID,
	reference = 'mitoLoci_' + fastafile + '.fasta', #your reference or reference sequences
	#PARAMETERS for NOVOALIGN
	insertSize='160', #mean fragment length of your libraries
	insertSTDEV = '60', #standard deviation of length of your libraries
	aScore = '90', #stringency score...I don't remember what 90 means, but it is stringent
	a1 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC', #adapter sequence
	a2 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA', #adapter sequence
	truncatePaired ='100',
	truncateUnpaired ='200',
	read1 = ID + r1name,
	read2 = ID + r2name,
	unpaired = ID + uname,
	novoOutPaired = ID + '_out_paired',
	novoOutUnpaired = ID + '_out_unpaired',
	outfile = ID + '_mito_sorted') #name your output

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
	novoalign -d {reference}.ndx -f {read1} {read2} -i PE {insertSize} {insertSTDEV} -t {aScore} -n {truncatePaired} -a {a1} {a2} -F STDFQ -o SAM > {novoOutPaired}
	novoalign -d {reference}.ndx -f {unpaired} -t 120 -n {truncateUnpaired} -a {a1} {a2} -F STDFQ -o SAM > {novoOutUnpaired}
	grep -v ZS:Z:NM {novoOutPaired} > {novoOutPaired}.sam
	grep -v ZS:Z:NM {novoOutUnpaired} > {novoOutUnpaired}.sam
	samtools view -bS {novoOutPaired}.sam > {novoOutPaired}.bam
	samtools view -bS {novoOutUnpaired}.sam > {novoOutUnpaired}.bam
	samtools merge -f {sample}raw.bam {novoOutPaired}.bam {novoOutUnpaired}.bam
	samtools sort {sample}raw.bam {outfile}
	samtools index {outfile}.bam
	samtools mpileup -f {reference} {outfile}.bam > {outfile}.pileup
	""".format(**variables)

	cmd_list = commands.split("\n")
	for cmd in cmd_list:
		os.system(cmd)

def main():
	args = get_args() 

	#Make a list of lists, each list within the list will have the first and second elements of the map file that are separated by a tab
	mylist = []
	with open(args.map) as rfile:
		for line in rfile:
			line = line.strip()
			mylist.append(line.split("\t"))

	#start the multiprocessing
	pool = multiprocessing.Pool()
	pool.map(align, mylist)#run the function with the arguments

if __name__ == "__main__": #run main over multiple processors
	main()







