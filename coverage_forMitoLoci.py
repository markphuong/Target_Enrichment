#!/usr/bin/env python

#SCRIPT TO: pull mtDNA genes from mtDNA genome

from __future__ import division
from collections import defaultdict
import os
import sys
import argparse
import multiprocessing
import commands


def get_args(): #arguments needed to give to this script
	parser = argparse.ArgumentParser(description="run coverage")

	#forces required argument to let it run
	required = parser.add_argument_group("required arguments") 
	required.add_argument("--map", help="textfile with samples to run and what fasta file to match it to", required=True) #A map file with the sample ID and the fasta file it goes to
	required.add_argument("--loci", help="textfile with desired loci", required=True)
	return parser.parse_args()

def coverage(element):

	index = element[0]
	lineage = element[1]


	pileupfile = index + '_mito_sorted.pileup'
	blastoutputfile = locus + '_mitoBlastOutput'
	ID = index + '_' + locus + '_' + lineage
	outfile = index + '_' + locus

	pileup = open(pileupfile, 'r')
	mypileup = defaultdict(list)
	for line in pileup:
		line = line.strip()
		line = line.split("\t")
		mypileup[line[0]].append(line[3]) #gene -> depth
	
	
	blastFile = open(blastoutputfile, 'r')
	for line in blastFile:
		info = line.strip()
		info = info.split("\t")
		if info[1] == locus and \
			info[0] == ID: 

			start = int(info[6]) - 1
			end = int(info[7]) 
			total = end-start
			total = float(total)
			numbers = mypileup[locus][start:end]
			#print numbers
			numbers = [float(i) for i in numbers]
			answer = sum(numbers)
			finalAnswer = answer/total
			
			outanswer = str(finalAnswer) + '\n'
			outSample = index + '_' + locus + '\t'
			out = open(outfile, 'w')
			out.write (outSample)
			out.write(outanswer)
			out.close
			break

def main():
	args = get_args() 

	#Make a list of lists, each list within the list will have the first and second elements of the map file that are separated by a tab
	mylist = []
	with open(args.map) as rfile:
		for line in rfile:
			line = line.strip()
			mylist.append(line.split("\t"))

	pool = multiprocessing.Pool()
	pool.map(coverage, mylist)#run the function with the arguments

args = get_args()
with open(args.loci) as rfile:
	for line in rfile:
		locus = line.strip()
		if __name__ == "__main__": #run main over multiple processors
			main()










