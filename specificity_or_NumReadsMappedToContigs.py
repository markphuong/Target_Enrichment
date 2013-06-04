#!/usr/bin/env python

#calculate specificity or, number of reads mapped to targets

from __future__ import division
import os
import sys
import argparse
import multiprocessing
import commands

def get_args(): #arguments needed to give to this script
	parser = argparse.ArgumentParser(description="run specificity")

	#forces required argument to let it run
	required = parser.add_argument_group("required arguments") 
	required.add_argument("--map", help="textfile with samples to run and what fasta file to match it to", required=True) #A map file with the sample ID and the fasta file it goes to

	return parser.parse_args()

def specificity(element):

	ID = element
	
	r1name = '_1p_final_renamed.fastq'  #extension of front reads
	r2name = '_2p_final_renamed.fastq' #extension of back reads
	uname = '_u_final_u_combined.fastq' #extension of unpaired reads

	variables = dict(
	sample = ID,
	read1 = ID + r1name,
	read2 = ID + r2name,
	unpaired = ID + uname,
	outPaired = ID + '_out_paired.sam',
	outUnpaired = ID + '_out_unpaired.sam') 

	command = """grep -c "HS" {outPaired}
grep -c "HS" {outUnpaired}
grep -c "@HS1" {read1}
grep -c "@HS1" {read2}
grep -c "@HS1" {unpaired}""".format(**variables)
	
	mynum = []	

	cmd_list = command.split("\n")
	for cmd in cmd_list:
		status, output = commands.getstatusoutput(cmd)
		mynum.append(float(output))

	answer = (mynum[0]+mynum[1]) / (mynum[2] + mynum[3] + mynum[4])*100

	outname = ID + '\t'
	outanswer = str(answer) + '\n'
	out.write(outname)
	out.write(outanswer)
	


	
args = get_args() 
outfile = 'specificity.txt'
out = open(outfile, 'w')
#Make a list of lists, each list within the list will have the first and second elements of the map file that are separated by a tab
mylist = []
with open(args.map) as rfile:
	for line in rfile:
		line = line.strip()
		myinput = line.split("\t")
		print myinput[0]
		specificity(myinput[0])
	
out.close













